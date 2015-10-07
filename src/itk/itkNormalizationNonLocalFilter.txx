/* ----------------------------- MNI Header -----------------------------------
@COPYRIGHT  :
              Copyright 2014 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

              This program is free software: you can redistribute it and/or modify
              it under the terms of the GNU General Public License as published by
              the Free Software Foundation, either version 3 of the License, or
              (at your option) any later version.
---------------------------------------------------------------------------- */
#ifndef __mincNormalizationNonLocalFilter_txx
#define __mincNormalizationNonLocalFilter_txx

namespace itk {
  template<class TFeatureImage, class TOutputImage, class TSearch, class TPatch,class SampleProcess,class TDistance,class TPreselectionFilter,class TRealType>
    std::ostream & operator<<(std::ostream &os, const typename NormalizationNonLocalFilter<TFeatureImage, TOutputImage, TSearch, TPatch,SampleProcess, TDistance,TPreselectionFilter,TRealType>::NormalizationLibraryType &s)
    {
      os<< " [ NormalizationLibraryType , size="<<s.size()<<" ]";
    }
    
  template<class TFeatureImage, class TOutputImage, class TSearch, class TPatch,class SampleProcess,class TDistance,class TPreselectionFilter,class TRealType>
  NormalizationNonLocalFilter<TFeatureImage, TOutputImage, TSearch, TPatch, SampleProcess,TDistance,TPreselectionFilter,TRealType>
  ::NormalizationNonLocalFilter():
    m_SearchKernel(),m_PatchKernel()
  {
    m_DefaultBoundaryCondition.SetConstant( itk::NumericTraits<PixelType>::Zero );
    m_BoundaryCondition = &m_DefaultBoundaryCondition;
  }

  template<class TFeatureImage, class TOutputImage, class TSearch,class TPatch,class SampleProcess,class TDistance,class TPreselectionFilter,class TRealType>
  void
  NormalizationNonLocalFilter<TFeatureImage, TOutputImage, TSearch, TPatch,SampleProcess, TDistance,TPreselectionFilter,TRealType>
  ::GenerateInputRequestedRegion()
  {
    // call the superclass' implementation of this method
    Superclass::GenerateInputRequestedRegion();
    
    // get pointers to the input and output
    typename Superclass::InputImagePointer  inputPtr = 
      const_cast< TFeatureImage * >( this->GetInput() );
    
    if ( !inputPtr )
    {
      return;
    }

    // get a copy of the input requested region (should equal the output
    // requested region)

    typename TFeatureImage::RegionType inputRequestedRegion;
    inputRequestedRegion = inputPtr->GetRequestedRegion();

    // pad the input requested region by the operator radius
    inputRequestedRegion.PadByRadius( m_SearchKernel.GetRadius() + m_PatchKernel.GetRadius() );

    // crop the input requested region at the input's largest possible region
    if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()) )
    {
      inputPtr->SetRequestedRegion( inputRequestedRegion );
      return;
    } else {
      // Couldn't crop the region (requested region is outside the largest
      // possible region).  Throw an exception.

      // store what we tried to request (prior to trying to crop)
      inputPtr->SetRequestedRegion( inputRequestedRegion );
      
      // build an exception
      itk::InvalidRequestedRegionError e(__FILE__, __LINE__);
      e.SetLocation(ITK_LOCATION);
      e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
      e.SetDataObject(inputPtr);
      throw e;
    }
  }

  template<class TFeatureImage, class TOutputImage, class TSearch,class TPatch,class SampleProcess,class TDistance,class TPreselectionFilter,class TRealType>
  void 
  NormalizationNonLocalFilter<TFeatureImage, TOutputImage, TSearch, TPatch,SampleProcess, TDistance,TPreselectionFilter,TRealType>
  ::BeforeThreadedGenerateData(void)
  {
    //check if the segmentation library is not empty, and that all the images have same dimensions
    if(m_NormalizationLibrary.empty())
    {
      itk::DataObjectError e(__FILE__,__LINE__);
      e.SetDescription("Normalization library is empty!");
      throw e;
    }

    typename Superclass::InputImagePointer  inputPtr = const_cast< TFeatureImage * >( this->GetInput() );

    for(size_t i=0;i<m_NormalizationLibrary.size();i++)
    {
      if( m_NormalizationLibrary[i]->GetLargestPossibleRegion()!=inputPtr->GetLargestPossibleRegion() )
      {
        itk::DataObjectError e(__FILE__,__LINE__);
        e.SetDescription("Normalization library image size mismatch!");
        throw e;
      }
    }

    if(m_PreloadedResults.IsNotNull())
    {
      this->GraftOutput(m_PreloadedResults);
    }
  }

  template<class TFeatureImage, class TOutputImage, class TSearch,
          class TPatch,class SampleProcess,class TDistance,class TPreselectionFilter,class TRealType>
  void
  NormalizationNonLocalFilter<TFeatureImage, TOutputImage, TSearch, 
                            TPatch,SampleProcess, TDistance,TPreselectionFilter,TRealType>
  ::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                        itk::ThreadIdType threadId) 
  {
    // Neighborhood iterators
    NeighborhoodIteratorType search_iter;

    //patch Neighborhood iterators
    NeighborhoodIteratorType patch_iter;

    // Find the boundary "faces"
    typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<FeatureImageType>::FaceListType faceList;

    itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<FeatureImageType> fC;

    faceList = fC(this->GetInput(0), outputRegionForThread, m_SearchKernel.GetRadius());

    typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<FeatureImageType>::FaceListType::iterator fit;

    itk::ImageRegionIterator<TOutputImage> o_iter;
    itk::ImageRegionConstIterator<RoiImageType> roi_iterator;

    itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

    // Process the boundary faces, these are N-d regions which border the
    // edge of the buffer

    const SearchKernelIteratorType searchKernelBegin = m_SearchKernel.Begin();
    const SearchKernelIteratorType searchKernelEnd   = m_SearchKernel.End();

    const PatchKernelIteratorType patchKernelBegin = m_PatchKernel.Begin();
    const PatchKernelIteratorType patchKernelEnd   = m_PatchKernel.End();

    //TODO: make a clone here ?
    PreselectionFilterPointerType preselection_filter=m_PreselectionFilter; 

    size_t library_size=m_NormalizationLibrary.size();

    std::vector<NeighborhoodIteratorType> patch_iterators(library_size);

    NormalizationSampleDistance sample_distance;

    sample_distance.reserve(library_size);

    for (fit = faceList.begin() ; fit != faceList.end(); ++fit)
    {
      search_iter = NeighborhoodIteratorType(m_SearchKernel.GetRadius(),
                                        this->GetInput(0), *fit);

      patch_iter = NeighborhoodIteratorType(m_PatchKernel.GetRadius(),
                                        this->GetInput(0), *fit);

      for(size_t sample=0;sample<library_size;sample++)
      {
        patch_iterators[sample] = NeighborhoodIteratorType(m_PatchKernel.GetRadius(),
                                        m_NormalizationLibrary[sample], *fit);

        patch_iterators[sample].OverrideBoundaryCondition(m_BoundaryCondition);
      }

      // primary output
      o_iter    = itk::ImageRegionIterator<OutputImageType>(this->GetOutput(), *fit);
      patch_iter.OverrideBoundaryCondition(m_BoundaryCondition);
      o_iter.GoToBegin();
      
      if(m_RoiImage.IsNotNull())
      {
        roi_iterator = itk::ImageRegionConstIterator<RoiImageType>(m_RoiImage,*fit);
        roi_iterator.GoToBegin();
      }
      
      
      while ( ! o_iter.IsAtEnd() )
      {
        bool within_roi=true;
        
        if(m_RoiImage.IsNotNull())
        {
          if(roi_iterator.Get()==0)
          {
            within_roi=false;
            o_iter.Set(m_Process->default_value());
          }
          ++roi_iterator;
        }

        if(!within_roi ) 
        {
          ++patch_iter;
          ++search_iter;
          
          ++o_iter;

          progress.CompletedPixel();
          continue;
        }

        sample_distance.clear();

        //we want to search closest-matching value from each sample
        for(size_t sample=0;sample<library_size;sample++)
        {
          double smallest_distance=-1.0;
          double best_sample=0.0;

          SearchKernelIteratorType kernel_it;
          ::size_t i;
          for(i=0, kernel_it=searchKernelBegin; kernel_it<searchKernelEnd; ++kernel_it, ++i )
          {

            // if structuring element is positive, use the pixel under that element
            // in the image
            if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero )
            {
              patch_iterators[ sample ].SetLocation(search_iter.GetIndex(i)); //move patch

              if(!patch_iterators[ sample ].InBounds()) 
                continue;

              if(!preselection_filter->select( patch_iter.GetIndex(),sample,patch_iterators[sample].GetIndex())) 
                continue;

              double dist=m_Distance->distance(patch_iter, patch_iterators[sample], patchKernelBegin, patchKernelEnd);

              if(smallest_distance<0.0 || dist<smallest_distance)
              {
                smallest_distance=dist;
                best_sample=::log(fabs(patch_iterators[sample].GetPixel(i)));
              }
            }

          }

          sample_distance.push_back(  NormSampleDistance<PixelType>(best_sample,smallest_distance,sample ) );

        }

        OutputPixelType output=0.0;

        m_Process->process(sample_distance, output );

        //TODO: maybe have a templated functor here too?
        o_iter.Set( ::log( fabs( patch_iter.GetPixel(patch_iter.Size() / 2) ) ) - output );
        //o_iter.Set( output );

        ++search_iter;
        ++o_iter;
        ++patch_iter;

        progress.CompletedPixel();

      }
    }
  }


  template<class TFeatureImage, class TOutputImage, class TSearch,class TPatch,class SampleProcess,class TDistance,class TPreselectionFilter,class TRealType>
  void 
  NormalizationNonLocalFilter<TFeatureImage, TOutputImage, TSearch, TPatch,SampleProcess, TDistance,TPreselectionFilter,TRealType>
  ::AfterThreadedGenerateData(void)
  {
  }

  template<class TFeatureImage, class TOutputImage, class TSearch,class TPatch,class SampleProcess,class TDistance,class TPreselectionFilter,class TRealType>
  void NormalizationNonLocalFilter<TFeatureImage, TOutputImage, TSearch, TPatch,SampleProcess, TDistance,TPreselectionFilter,TRealType>
  ::PrintSelf(std::ostream &os, itk::Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
    
    os << indent << "Search Kernel: " << m_SearchKernel << std::endl;
    os << indent << "Patch Kernel: "  << m_PatchKernel  << std::endl;
    os << indent << "Boundary condition: " << typeid( *m_BoundaryCondition ).name() << std::endl;
    
  }

}// end namespace itk
#endif

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl C++
