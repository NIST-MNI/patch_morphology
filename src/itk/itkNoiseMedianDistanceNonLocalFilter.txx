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
#ifndef __mincNoiseMedianDistanceNonLocalFilter_txx
#define __mincNoiseMedianDistanceNonLocalFilter_txx
#include <itkProgressAccumulator.h>


namespace itk {

template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
NoiseMedianDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::NoiseMedianDistanceNonLocalFilter():
 m_SearchKernel(),m_PatchKernel()
{
  m_DefaultBoundaryCondition.SetConstant( itk::NumericTraits<PixelType>::Zero );
  m_BoundaryCondition = &m_DefaultBoundaryCondition;
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void 
NoiseMedianDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
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
NoiseMedianDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
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

    
  typename itk::SubtractImageFilter<TOutputImage,TOutputImage,TOutputImage>::Pointer 
    sub2 = itk::SubtractImageFilter<TOutputImage,TOutputImage,TOutputImage>::New();
    
  typename itk::DivideImageFilter<TOutputImage,TOutputImage,TOutputImage>::Pointer 
    div = itk::DivideImageFilter<TOutputImage,TOutputImage,TOutputImage>::New();
    
  typename MinimalDistanceNonLocalFilter<TInputImage,TOutputImage,TPatch,TSearch,TDistance>::Pointer 
    mdist = MinimalDistanceNonLocalFilter<TInputImage,TOutputImage,TPatch,TSearch,TDistance>::New();
    
  typename MedianDistanceNonLocalFilter<TInputImage,TOutputImage,TPatch,TSearch,TDistance>::Pointer 
    median_dist = MedianDistanceNonLocalFilter<TInputImage,TOutputImage,TPatch,TSearch,TDistance>::New();

    
 
  mean->SetInput( this->GetInput());
  mean->SetKernel( this->GetPatchKernel());
  
  mean->OverrideBoundaryCondition(m_BoundaryCondition);
  
  sub->SetInPlace(false);
  sub->SetInput1(this->GetInput());
  sub->SetInput2(mean->GetOutput());
  
  mdist->SetSearchKernel( this->GetSearchKernel() );
  mdist->SetPatchKernel( this->GetPatchKernel() );
  mdist->OverrideBoundaryCondition(m_BoundaryCondition);
  mdist->SetInput(sub->GetOutput());
  
  median_dist->SetSearchKernel( this->GetSearchKernel() );
  median_dist->SetPatchKernel( this->GetPatchKernel() );
  median_dist->OverrideBoundaryCondition(m_BoundaryCondition);
  median_dist->SetInput(this->GetInput());
  median_dist->SetPreselectionFilter(m_PreselectionFilter);
  
  div->SetInPlace(false);
  div->SetInput2(median_dist->GetOutput());
  div->SetInput1(mdist->GetOutput());
  div->GraftOutput( this->GetOutput() );
//   sub2->SetInPlace(false);
//   sub2->SetInput1(median_dist->GetOutput());
//   sub2->SetInput2(mdist->GetOutput());

  if(m_RoiImage.IsNotNull())
  {
    mdist->SetRoiImage(m_RoiImage);
    median_dist->SetRoiImage(m_RoiImage);
  }
  
  
  //sub2->GraftOutput( this->GetOutput() );
  
  double total=this->GetPatchKernel().Size()*this->GetSearchKernel().Size()*2+this->GetPatchKernel().Size()+1;
  
  progress->RegisterInternalFilter(mean,this->GetPatchKernel().Size()/total);
  progress->RegisterInternalFilter(sub,1.0/total);
  progress->RegisterInternalFilter(mdist,this->GetPatchKernel().Size()*this->GetSearchKernel().Size()/total);
  progress->RegisterInternalFilter(median_dist,this->GetPatchKernel().Size()*this->GetSearchKernel().Size()/total);
  progress->RegisterInternalFilter(div,1.0/total);
//  progress->RegisterInternalFilter(sub2,1.0/total);
  
  div->Update();
  //sub2->Update();
  this->GraftOutput( div->GetOutput() );
  //this->GraftOutput( sub2->GetOutput() );
}


template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void NoiseMedianDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  os << indent << "Search Kernel: " << m_SearchKernel << std::endl;
  os << indent << "Patch Kernel: "  << m_PatchKernel  << std::endl;
  os << indent << "Boundary condition: " << typeid( *m_BoundaryCondition ).name() << std::endl;
  os << indent << "Preselection filter: " << typeid( m_PreselectionFilter ).name() << std::endl;
  
}

}// end namespace itk
#endif

// kate: space-indent on; hl c++;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 2 