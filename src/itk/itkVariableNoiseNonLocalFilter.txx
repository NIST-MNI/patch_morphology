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
#ifndef __mincVariableNoiseNonLocalFilter_txx
#define __mincVariableNoiseNonLocalFilter_txx


namespace itk {

  
template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
VariableNoiseNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::VariableNoiseNonLocalFilter():
    m_SearchKernel(),
    m_PatchKernel(),
    m_OutputMeanWeight(false),
    m_Beta(1.0)
{
  m_DefaultBoundaryCondition.SetConstant( itk::NumericTraits<PixelType>::Zero );
  m_BoundaryCondition = &m_DefaultBoundaryCondition;
  this->SetNumberOfRequiredInputs(2);
  m_Distance = TDistance::New();
  m_PreselectionFilter = TPreselectionFilter::New();
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void 
VariableNoiseNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // get pointers to the input and output
  typename Superclass::InputImagePointer  inputPtr = 
    const_cast< TInputImage * >( this->GetInput(0) );
  
  if ( !inputPtr )
  {
    return;
  }
  
  // get pointers to the input and output
  typename Superclass::InputImagePointer  inputPtr2 = 
    const_cast< TInputImage * >( this->GetInput(1) );
  
  if ( !inputPtr2 )
  {
    return;
  }

  // get a copy of the input requested region (should equal the output
  // requested region)

  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius( m_SearchKernel.GetRadius() + m_PatchKernel.GetRadius() );

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()) )
  {
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    inputPtr2->SetRequestedRegion( inputRequestedRegion );
    return;
  } else {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    inputPtr2->SetRequestedRegion( inputRequestedRegion );
    
    // build an exception
    itk::InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
  }
  
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void 
VariableNoiseNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::BeforeThreadedGenerateData(void)
{
  this->m_Weight.SetBeta(m_Beta);
}



#if (ITK_VERSION_MAJOR==3)
template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void
VariableNoiseNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                      int threadId) 
#else
template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void
VariableNoiseNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                      itk::ThreadIdType threadId) 
#endif    

{
  // Neighborhood iterators
  NeighborhoodIteratorType search_iter;
  
  //patch Neighborhood iterators
  NeighborhoodIteratorType patch_iter1;
  NeighborhoodIteratorType patch_iter2;

  // Find the boundary "faces"
  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList,faceList2;
  
  itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> fC;
  itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> fC2;
  
  faceList = fC(this->GetInput(0), outputRegionForThread, m_SearchKernel.GetRadius());
  faceList2 = fC(this->GetInput(1), outputRegionForThread, m_SearchKernel.GetRadius());

  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit1,fit2;

  itk::ImageRegionIterator<TOutputImage> o_iter;
  itk::ImageRegionConstIterator<TInputImage> i_sigma2;
  itk::ImageRegionConstIterator<RoiImageType> roi_iterator;

  itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  // Process the boundary faces, these are N-d regions which border the
  // edge of the buffer

  const SearchKernelIteratorType searchKernelBegin = m_SearchKernel.Begin();
  const SearchKernelIteratorType searchKernelEnd   = m_SearchKernel.End();

  const SearchKernelIteratorType patchKernelBegin = m_PatchKernel.Begin();
  const SearchKernelIteratorType patchKernelEnd   = m_PatchKernel.End();
  PreselectionFilterPointerType preselection=m_PreselectionFilter;

  for (fit1 = faceList.begin(),fit2 = faceList2.begin() ; fit1 != faceList.end(); ++fit1, ++fit2)
  { 
    search_iter = NeighborhoodIteratorType(m_SearchKernel.GetRadius(),
                                      this->GetInput(0), *fit1);
    
    patch_iter1 = NeighborhoodIteratorType(m_PatchKernel.GetRadius(),
                                      this->GetInput(0), *fit1);

    patch_iter2 = NeighborhoodIteratorType(m_PatchKernel.GetRadius(),
                                      this->GetInput(0), *fit1);
    
    o_iter = itk::ImageRegionIterator<OutputImageType>(this->GetOutput(), *fit1);
    i_sigma2= itk::ImageRegionConstIterator<OutputImageType>(this->GetInput(1), *fit2);
    
    search_iter.OverrideBoundaryCondition(m_BoundaryCondition);
    patch_iter1.OverrideBoundaryCondition(m_BoundaryCondition);
    patch_iter2.OverrideBoundaryCondition(m_BoundaryCondition);

    search_iter.GoToBegin();
    o_iter.GoToBegin();
    i_sigma2.GoToBegin();
    
    if(m_RoiImage.IsNotNull())
    {
      roi_iterator = itk::ImageRegionConstIterator<RoiImageType>(m_RoiImage,*fit1);
      roi_iterator.GoToBegin();
      
      while ( ! o_iter.IsAtEnd() )
      {
        if(roi_iterator.Get()==0)
        {
          ++roi_iterator;
          ++search_iter;
          o_iter.Set(0); //TODO: maybe put input value here?
          ++o_iter;
          ++patch_iter1;
          progress.CompletedPixel();
          continue;
        } 
        ++roi_iterator;
        
        o_iter.Set( this->Evaluate( search_iter, (double)i_sigma2.Get(), patch_iter1, patch_iter2, searchKernelBegin, searchKernelEnd, patchKernelBegin, patchKernelEnd,preselection) );
        ++search_iter;
        ++o_iter;
        ++patch_iter1;
        ++i_sigma2;
        progress.CompletedPixel();
      }
    }	else {	
      while ( ! o_iter.IsAtEnd() )
      {
        o_iter.Set( this->Evaluate( search_iter, (double)i_sigma2.Get(), patch_iter1, patch_iter2, searchKernelBegin, searchKernelEnd, patchKernelBegin, patchKernelEnd,preselection) );
        ++search_iter;
        ++o_iter;
        ++patch_iter1;
        ++i_sigma2;
        progress.CompletedPixel();
      }
    }

  }
}


template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
typename VariableNoiseNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>::PixelType
VariableNoiseNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
  ::Evaluate(const NeighborhoodIteratorType &searchIt,
              double sigma2,
              const NeighborhoodIteratorType &patchIt1,
              NeighborhoodIteratorType &patchIt2,
              const SearchKernelIteratorType searchKernelBegin,
              const SearchKernelIteratorType searchKernelEnd,
              const PatchKernelIteratorType patchKernelBegin,
              const PatchKernelIteratorType patchKernelEnd,
              PreselectionFilterPointerType flt
            )
{
  static const double SIGMA_EPSILON=1e-6;  
  
  ::size_t i;
  SearchKernelIteratorType kernel_it;
  ::size_t center = (::size_t) (searchIt.Size() / 2); // get offset of center pixel

  double total_weight=0.0;
  double total=0.0;
  double smallest_weight=0.0;
  //variable noise part 
  
  if( sigma2 > SIGMA_EPSILON )
  {
    for( i=0, kernel_it=searchKernelBegin; kernel_it<searchKernelEnd; ++kernel_it, ++i )
    {
      // if structuring element is positive, use the pixel under that element
      // in the image
      if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero )
      {
        patchIt2.SetLocation(searchIt.GetIndex(i)); //move patch
        
        if(!patchIt2.InBounds()) continue;
        if(!flt->select(patchIt1.GetIndex(),patchIt2.GetIndex())) continue;
        
        double distance=0.0;
        
        if(i!=center) 
          distance=m_Distance->distance(patchIt1,patchIt2,patchKernelBegin,patchKernelEnd);

        double weight=m_Weight(distance,sigma2);
        
        total_weight+=weight;
        total+=searchIt.GetPixel(i)*weight;
      }
    }
  }
  
  if(total_weight==0.0)
  {
    total=searchIt.GetPixel(center);
    total_weight=1.0;
  }
  
  if(m_OutputMeanWeight)
    return (PixelType)(total_weight/i);
  else
    return (PixelType)(total/total_weight);
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void VariableNoiseNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  os << indent << "Search Kernel: "      << m_SearchKernel << std::endl;
  os << indent << "Patch Kernel: "       << m_PatchKernel  << std::endl;
  os << indent << "Boundary condition: " << typeid( *m_BoundaryCondition ).name() << std::endl;
  os << indent << "Beta:"                << m_Beta << std::endl;
  
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void VariableNoiseNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::SetInput1( const TInputImage * image1)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<TInputImage *>( image1 ));
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void VariableNoiseNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::SetInput2( const TInputImage * image2)
{
  this->SetNthInput(1, const_cast<TInputImage *>( image2 ));
}


}// end namespace itk
#endif

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl C++
