#ifndef __mincRegressionNonLocalFilter_txx
#define __mincRegressionNonLocalFilter_txx

#include <set>
 
namespace itk {

  
template<class TFeatureImage,class TLabelImage, class TOutputImage, class TSearch, class TPatch,class TRegressor,class TPreselectionFilter,class TRealType>
  std::ostream & operator<<(std::ostream &os, const typename RegressionNonLocalFilter<TFeatureImage, TLabelImage, TOutputImage, TSearch, TRegressor, TPreselectionFilter, TRealType>::SegmentationLibraryType &s)
  {
    os<< " [ RegressionLibraryType , size="<<s.size()<<" ]";
  }
  
template<class TFeatureImage,class TLabelImage, class TOutputImage, class TSearch, class TPatch,class TRegressor,class TPreselectionFilter,class TRealType>
RegressionNonLocalFilter<TFeatureImage,TLabelImage,TOutputImage,TSearch,TPatch,TRegressor,TPreselectionFilter,TRealType>
::RegressionNonLocalFilter():
  m_SearchKernel(),
  m_PatchKernel(),
  m_LabelCount(2)
{
  m_BoundaryCondition = &m_DefaultBoundaryCondition;
  
  m_DefaultLabelBoundaryCondition.SetConstant( itk::NumericTraits<LabelPixelType>::Zero );
  m_LabelBoundaryCondition=&m_DefaultLabelBoundaryCondition;
  
  m_ConfidenceThreshold=0.5; //Run everywhere where confidence <0.5
  m_NonConfidentVoxels=0;
}

template<class TFeatureImage,class TLabelImage, class TOutputImage, class TSearch, class TPatch,class TRegressor,class TPreselectionFilter,class TRealType>
void
RegressionNonLocalFilter<TFeatureImage,TLabelImage,TOutputImage,TSearch,TPatch,TRegressor,TPreselectionFilter,TRealType>
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

template<class TFeatureImage,class TLabelImage, class TOutputImage, class TSearch, class TPatch,class TRegressor,class TPreselectionFilter,class TRealType>
void 
RegressionNonLocalFilter<TFeatureImage,TLabelImage,TOutputImage,TSearch,TPatch,TRegressor,TPreselectionFilter,TRealType>
::BeforeThreadedGenerateData(void)
{
  //check if the segmentation library is not empty, and that all the images have same dimensions
  if(m_SegmentationLibrary.empty())
  {
    itk::DataObjectError e(__FILE__,__LINE__);
    e.SetDescription("Segmentation library is empty!");
    throw e;
  }

  typename Superclass::InputImagePointer  inputPtr = const_cast< TFeatureImage * >( this->GetInput() );

  for(size_t i=0;i<m_SegmentationLibrary.size();i++)
  {
    if( m_SegmentationLibrary[i].feature->GetLargestPossibleRegion()!=inputPtr->GetLargestPossibleRegion() ||
        m_SegmentationLibrary[i].label->GetLargestPossibleRegion()!=inputPtr->GetLargestPossibleRegion() )
    {
      itk::DataObjectError e(__FILE__,__LINE__);
      e.SetDescription("Segmentation library image size mismatch!");
      throw e;
    }
  }

  if(m_Confidence.IsNull())
  {
    m_Confidence=FloatImageType::New();
    //allocate space in the confidence field and search field
    m_Confidence->SetLargestPossibleRegion(inputPtr->GetLargestPossibleRegion());
    m_Confidence->SetBufferedRegion(inputPtr->GetLargestPossibleRegion());
    m_Confidence->SetRequestedRegion(inputPtr->GetLargestPossibleRegion());
    m_Confidence->SetSpacing( inputPtr->GetSpacing() );
    m_Confidence->SetOrigin ( inputPtr->GetOrigin() );
    m_Confidence->SetDirection(inputPtr->GetDirection());
    m_Confidence->Allocate();
    m_Confidence->FillBuffer(-1.0);
  } else { //reuse old confidence results
    if(m_Confidence->GetLargestPossibleRegion()!=inputPtr->GetLargestPossibleRegion())
    {
      itk::DataObjectError e(__FILE__,__LINE__);
      e.SetDescription("Confidence image size mismatch!");
      throw e;
    }
  }

  if(m_SearchDistance.IsNull())
  {
    m_SearchDistance=FloatImageType::New();
    m_SearchDistance->SetLargestPossibleRegion(inputPtr->GetLargestPossibleRegion());
    m_SearchDistance->SetBufferedRegion(inputPtr->GetLargestPossibleRegion());
    m_SearchDistance->SetRequestedRegion(inputPtr->GetLargestPossibleRegion());
    m_SearchDistance->SetSpacing( inputPtr->GetSpacing() );
    m_SearchDistance->SetOrigin ( inputPtr->GetOrigin() );
    m_SearchDistance->SetDirection(inputPtr->GetDirection());
    m_SearchDistance->Allocate();
    m_SearchDistance->FillBuffer(0.0);
  } else { //reuse old confidence results
    if(m_SearchDistance->GetLargestPossibleRegion()!=inputPtr->GetLargestPossibleRegion())
    {
      itk::DataObjectError e(__FILE__,__LINE__);
      e.SetDescription("SearchDistance image size mismatch!");
      throw e;
    }
  }
  
  if(m_GradingMap.IsNull())
  {
    m_GradingMap=FloatImageType::New();
    m_GradingMap->SetLargestPossibleRegion(inputPtr->GetLargestPossibleRegion());
    m_GradingMap->SetBufferedRegion(inputPtr->GetLargestPossibleRegion());
    m_GradingMap->SetRequestedRegion(inputPtr->GetLargestPossibleRegion());
    m_GradingMap->SetSpacing( inputPtr->GetSpacing() );
    m_GradingMap->SetOrigin ( inputPtr->GetOrigin() );
    m_GradingMap->SetDirection(inputPtr->GetDirection());
    m_GradingMap->Allocate();
    m_GradingMap->FillBuffer(0.0);
  } else { //reuse old confidence results
    if(m_GradingMap->GetLargestPossibleRegion()!=inputPtr->GetLargestPossibleRegion())
    {
      itk::DataObjectError e(__FILE__,__LINE__);
      e.SetDescription("SearchDistance image size mismatch!");
      throw e;
    }
  }
  
  if(m_ProbabilityMaps.empty())
  {
    m_ProbabilityMaps.resize(m_LabelCount);
    for(int i=0;i<m_LabelCount;i++)
    {
      if(m_ProbabilityMaps[i].IsNull())
      {
        m_ProbabilityMaps[i]=FloatImageType::New();
        m_ProbabilityMaps[i]->SetLargestPossibleRegion(inputPtr->GetLargestPossibleRegion());
        m_ProbabilityMaps[i]->SetBufferedRegion(inputPtr->GetLargestPossibleRegion());
        m_ProbabilityMaps[i]->SetRequestedRegion(inputPtr->GetLargestPossibleRegion());
        m_ProbabilityMaps[i]->SetSpacing( inputPtr->GetSpacing() );
        m_ProbabilityMaps[i]->SetOrigin ( inputPtr->GetOrigin() );
        m_ProbabilityMaps[i]->SetDirection(inputPtr->GetDirection());
        m_ProbabilityMaps[i]->Allocate();
        m_ProbabilityMaps[i]->FillBuffer(0.0);
      }
    }
  } else {
    if (m_ProbabilityMaps.size()!=m_LabelCount)
    {
      itk::DataObjectError e(__FILE__,__LINE__);
      e.SetDescription("Wrong number of probability maps!");
      throw e;
    }
    for(int i=0;i<m_LabelCount;i++)
    {
      if(m_ProbabilityMaps[i]->GetLargestPossibleRegion()!=inputPtr->GetLargestPossibleRegion())
      {
        itk::DataObjectError e(__FILE__,__LINE__);
        e.SetDescription("probability image size mismatch!");
        throw e;
      }
    }
  }
  
  // Create the thread temporaries
  m_ThreadNonConfidentVoxels = std::vector<TRealType>(this->GetNumberOfThreads(),0);
  
  if(m_PreloadedResults.IsNotNull())
  {
    this->GraftOutput(m_PreloadedResults);
  }
}

template<class TFeatureImage,class TLabelImage, class TOutputImage, class TSearch, class TPatch,class TRegressor,class TPreselectionFilter,class TRealType>
void
RegressionNonLocalFilter<TFeatureImage,TLabelImage,TOutputImage,TSearch,TPatch,TRegressor,TPreselectionFilter,TRealType>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                      itk::ThreadIdType threadId)
{
  // disable dynamic multithreading
  this->DynamicMultiThreadingOff();
  
  // clone regressor , because it's not thread safe
  RegressorPointer _regressor=m_Regressor->clone();

  // Neighborhood iterators
  NeighborhoodIteratorType search_iter;
  
  //patch Neighborhood iterators
  NeighborhoodIteratorType patch_iter;
  
  //probability vector
  std::vector<double> prob(m_LabelCount,0.0);
  std::set<int> labels_present;

  // Find the boundary "faces"
  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<FeatureImageType>::FaceListType faceList;
  
  itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<FeatureImageType> fC;
  
  faceList = fC(this->GetInput(0), outputRegionForThread, m_SearchKernel.GetRadius());

  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<FeatureImageType>::FaceListType::iterator fit;

  itk::ImageRegionIterator<TOutputImage> o_iter;
  itk::ImageRegionConstIterator<RoiImageType> roi_iterator;
  
  itk::ImageRegionIterator<FloatImageType> conf_iter;
  itk::ImageRegionIterator<FloatImageType> sdist_iter;
  itk::ImageRegionIterator<FloatImageType> grading_iter;
  std::vector<itk::ImageRegionIterator<FloatImageType> > prob_iter(m_LabelCount);

  itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  // Process the boundary faces, these are N-d regions which border the
  // edge of the buffer

  const SearchKernelIteratorType searchKernelBegin = m_SearchKernel.Begin();
  const SearchKernelIteratorType searchKernelEnd   = m_SearchKernel.End();

  const SearchKernelIteratorType patchKernelBegin = m_PatchKernel.Begin();
  const SearchKernelIteratorType patchKernelEnd   = m_PatchKernel.End();
  
  //TODO: make a clone here ?
  PreselectionFilterPointerType preselection_filter=m_PreselectionFilter; 
  
  size_t library_size=m_SegmentationLibrary.size();
  

  std::vector<NeighborhoodIteratorType> patch_iterators(library_size);
  std::vector<LabelNeighborhoodIteratorType> search_iterators(library_size);
  
  for (fit = faceList.begin() ; fit != faceList.end(); ++fit)
  { 
    search_iter = NeighborhoodIteratorType(m_SearchKernel.GetRadius(),
                                      this->GetInput(0), *fit);
    
    patch_iter = NeighborhoodIteratorType(m_PatchKernel.GetRadius(),
                                      this->GetInput(0), *fit);
    
  
    for(size_t sample=0;sample<library_size;sample++)
    {
      patch_iterators[sample] = NeighborhoodIteratorType(m_PatchKernel.GetRadius(),
                                      m_SegmentationLibrary[sample].feature, *fit);
      
      patch_iterators[sample].OverrideBoundaryCondition(m_BoundaryCondition);
      
      search_iterators[sample] = LabelNeighborhoodIteratorType(m_SearchKernel.GetRadius(),
                                      m_SegmentationLibrary[sample].label, *fit);
      
      search_iterators[sample].OverrideBoundaryCondition(m_LabelBoundaryCondition);
      search_iterators[sample].GoToBegin();
    }
    
    // primary output
    o_iter    = itk::ImageRegionIterator<OutputImageType>(this->GetOutput(), *fit);
    
    //confidence interval output
    conf_iter = itk::ImageRegionIterator<FloatImageType>(m_Confidence, *fit);
    
    //mean search distance output
    sdist_iter= itk::ImageRegionIterator<FloatImageType>(m_SearchDistance, *fit);
    
    //grading output
    grading_iter= itk::ImageRegionIterator<FloatImageType>(m_GradingMap, *fit);

    search_iter.OverrideBoundaryCondition(m_BoundaryCondition);
    patch_iter.OverrideBoundaryCondition(m_BoundaryCondition);

    search_iter.GoToBegin();
    o_iter.GoToBegin();
    conf_iter.GoToBegin();
    sdist_iter.GoToBegin();
    grading_iter.GoToBegin();
    
    for(int l=0;l<m_LabelCount;l++)
    {
      prob_iter[l]=itk::ImageRegionIterator<FloatImageType>(m_ProbabilityMaps[l], *fit);
      prob_iter[l].GoToBegin();
    }
    
    if(m_RoiImage.IsNotNull())
    {
      roi_iterator = itk::ImageRegionConstIterator<RoiImageType>(m_RoiImage,*fit);
      roi_iterator.GoToBegin();
    }
    
    _regressor->allocate_memory(search_iter.Size()*library_size, patch_iter.Size(), m_LabelCount);

    while ( ! o_iter.IsAtEnd() )
    {
      bool within_roi=true;
      
      if(m_RoiImage.IsNotNull())
      {
        if(roi_iterator.Get()==0)
        {
          within_roi=false;
          o_iter.Set(_regressor->default_value());
        }
        ++roi_iterator;
      }
      
      if(!within_roi || conf_iter.Get() >= m_ConfidenceThreshold ) 
      {
        ++search_iter;
        ++o_iter;
        ++patch_iter;
        ++conf_iter;
        ++sdist_iter;
        ++grading_iter;
        
        for(int l=0;l<m_LabelCount;l++)
          ++prob_iter[l];
        
        for(size_t sample=0;sample<library_size;sample++)
          ++search_iterators[sample]; 

        progress.CompletedPixel();
        continue;
      }
      
      _regressor->start_new_regression();
      _regressor->set_input(patch_iter,patchKernelBegin, patchKernelEnd);
      labels_present.clear();
      //Iterate through all library samples and find distances
      for(size_t sample=0;sample<library_size;sample++)
      {
        SearchKernelIteratorType kernel_it;
        ::size_t center = (::size_t) (search_iter.Size() / 2); // get offset of center pixel

        double total_weight=0.0;
        double total=0.0;
        double smallest_weight=0.0;
        double sample_grading=m_SegmentationLibrary[sample].grading;
        
        ::size_t  i;
        
        
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
              
            labels_present.insert(search_iterators[sample].GetPixel(i));
            
            _regressor->add_training_sample(search_iterators[sample].GetPixel(i), // label or sample_grading
                                            sample_grading,
                                            patch_iterators[sample], //features
                                            patchKernelBegin, patchKernelEnd, // patch iterator
                                            sample //sample id
                                          );
          }
        }
        ++search_iterators[sample];
      }
      
      int    output=_regressor->default_value();
      double confidence=0;
      double search_distance=0;
      double grading=0;
      prob.assign(m_LabelCount,0.0);
      
      //TODO: make sure we are not trying to grade here!!!!!
      if(labels_present.size()<2)
      {
        output=*labels_present.begin();
        confidence=1;
        grading=-1;
        prob[output]=1.0;
      } else {
        if(!_regressor->regress(output,confidence,prob,grading))
        {
          //TODO: store information that regression failed on this voxel, somehow....
          output=0;
          confidence=0;
          
        }
      }
      //grading=labels_present.size();
      
      o_iter.Set( output );
      conf_iter.Set( confidence );
      sdist_iter.Set( search_distance );
      grading_iter.Set( grading );
      
      for(int l=0;l<m_LabelCount;l++)
      {
        prob_iter[l].Set( prob[l] );
        //prob_iter[l].Set( l );
        ++prob_iter[l];
      }
      
      ++search_iter;
      ++o_iter;
      ++patch_iter;
      ++conf_iter;
      ++sdist_iter;
      ++grading_iter;
      
      
      if(confidence<m_ConfidenceThreshold)
      {
        m_ThreadNonConfidentVoxels[threadId]+=1;
      }

      progress.CompletedPixel();
    }
  }
}


template<class TFeatureImage,class TLabelImage, class TOutputImage, class TSearch, class TPatch,class TRegressor,class TPreselectionFilter,class TRealType>
void 
RegressionNonLocalFilter<TFeatureImage,TLabelImage,TOutputImage,TSearch,TPatch,TRegressor,TPreselectionFilter,TRealType>
::AfterThreadedGenerateData(void)
{
  int numberOfThreads = this->GetNumberOfThreads();
  m_NonConfidentVoxels=0;
  for( int i = 0; i < numberOfThreads; i++)
    m_NonConfidentVoxels+=m_ThreadNonConfidentVoxels[i];
}

template<class TFeatureImage,class TLabelImage, class TOutputImage, class TSearch, class TPatch,class TRegressor,class TPreselectionFilter,class TRealType>
void
RegressionNonLocalFilter<TFeatureImage,TLabelImage,TOutputImage,TSearch,TPatch,TRegressor,TPreselectionFilter,TRealType>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  os << indent << "Search Kernel: " << m_SearchKernel << std::endl;
  os << indent << "Patch Kernel: "  << m_PatchKernel  << std::endl;
  os << indent << "Boundary condition: " << typeid( *m_BoundaryCondition ).name() << std::endl;
  
}

}// end namespace itk
#endif

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;show-tabs on;tab-width 2;hl c++
