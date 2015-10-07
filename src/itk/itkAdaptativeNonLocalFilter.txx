#ifndef __mincAdaptativeNonLocalFilter_txx
#define __mincAdaptativeNonLocalFilter_txx
#include <itkProgressAccumulator.h>


namespace itk {

template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
AdaptativeNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::AdaptativeNonLocalFilter():
 m_SearchKernel(),m_PatchKernel(),m_OutputMeanWeight(false),m_Beta(1.0),m_Regularize(0.0)
{
  m_DefaultBoundaryCondition.SetConstant( itk::NumericTraits<PixelType>::Zero );
  m_BoundaryCondition = &m_DefaultBoundaryCondition;
  m_Distance = TDistance::New();
  m_PreselectionFilter = TPreselectionFilter::New();
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void 
AdaptativeNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
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


template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void
AdaptativeNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::GenerateData() 
{
  // Create a process accumulator for tracking the progress of this minipipeline
  itk::ProgressAccumulator::Pointer progress = itk::ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  // Allocate the output
  this->AllocateOutputs();

  typename LocalMeanFilter<TInputImage,TOutputImage,TPatch>::Pointer 
    mean = LocalMeanFilter<TInputImage,TOutputImage,TPatch>::New();

  typename itk::SubtractImageFilter<TInputImage,TInputImage,TInputImage>::Pointer 
    sub = itk::SubtractImageFilter<TInputImage,TInputImage,TInputImage>::New();
    
  typename MinimalDistanceNonLocalFilter<TInputImage,TOutputImage,TPatch,TSearch,TDistance>::Pointer 
    mdist = MinimalDistanceNonLocalFilter<TInputImage,TOutputImage,TPatch,TSearch,TDistance>::New();
    
  typename VariableNoiseNonLocalFilter<TInputImage,TOutputImage,TPatch,TSearch,TDistance,TWeight,TPreselectionFilter>::Pointer 
    vn = VariableNoiseNonLocalFilter<TInputImage,TOutputImage,TPatch,TSearch,TDistance,TWeight,TPreselectionFilter>::New();

  typename itk::DiscreteGaussianImageFilter<TInputImage,TOutputImage>::Pointer 
    regularize = itk::DiscreteGaussianImageFilter<TInputImage,TOutputImage>::New();
    
  mean->SetInput( this->GetInput());
  mean->SetKernel( this->GetPatchKernel());

  mean->OverrideBoundaryCondition(m_BoundaryCondition);

  sub->SetInput1(this->GetInput());
  sub->SetInput2(mean->GetOutput());

  mdist->SetSearchKernel( this->GetSearchKernel() );
  mdist->SetPatchKernel( this->GetPatchKernel() );
  mdist->OverrideBoundaryCondition(m_BoundaryCondition);
  mdist->SetInput(sub->GetOutput());

  vn->SetSearchKernel( this->GetSearchKernel() );
  vn->SetPatchKernel( this->GetPatchKernel() );
  vn->SetInput1(this->GetInput());
  
  if(m_Regularize>0.1) 
  {
    regularize->SetInput(mdist->GetOutput());
    regularize->SetUseImageSpacing(true);
    regularize->SetVariance(m_Regularize);
    
    vn->SetInput2(regularize->GetOutput());
  } else {
    vn->SetInput2(mdist->GetOutput());
  }

  vn->OverrideBoundaryCondition(m_BoundaryCondition);
  vn->GraftOutput( this->GetOutput() );
  vn->SetPreselectionFilter(this->m_PreselectionFilter);
  vn->SetOutputMeanWeight(this->m_OutputMeanWeight);
  vn->SetBeta(this->m_Beta);

  double total=this->GetPatchKernel().Size()*this->GetSearchKernel().Size()*2.0
              +this->GetPatchKernel().Size()+1.0+(m_Regularize>0.1?1:0);

  progress->RegisterInternalFilter(mean,this->GetPatchKernel().Size()/total);
  progress->RegisterInternalFilter(sub,1.0/total);
  progress->RegisterInternalFilter(mdist,this->GetPatchKernel().Size()*this->GetSearchKernel().Size()/total);
  if(m_Regularize>0.1) 
  {
    progress->RegisterInternalFilter(regularize,1.0/total);
  }

  progress->RegisterInternalFilter(vn,this->GetPatchKernel().Size()*this->GetSearchKernel().Size()/total);

  if(m_RoiImage.IsNotNull())
  {
    mdist->SetRoiImage(m_RoiImage);
    vn->SetRoiImage(m_RoiImage);
  }

  vn->Update();
  this->GraftOutput( vn->GetOutput() );
}


template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void AdaptativeNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
	
  os << indent << "Search Kernel: " << m_SearchKernel << std::endl;
  os << indent << "Patch Kernel: "  << m_PatchKernel  << std::endl;
  os << indent << "Boundary condition: " << typeid( *m_BoundaryCondition ).name() << std::endl;
  os << indent << "Preselection filter: " << typeid( m_PreselectionFilter ).name() << std::endl;
  os << indent << "Beta: " << m_Beta << std::endl;
  os << indent << "Regularize: " << m_Regularize << std::endl;
	
}

}// end namespace itk
#endif

// kate: space-indent on; hl c++;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 2 
