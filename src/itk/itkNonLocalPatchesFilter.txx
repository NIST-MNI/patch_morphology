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
#ifndef __mincNonLocalPatchesFilter_txx
#define __mincNonLocalPatchesFilter_txx

#include <itkNeighborhoodAlgorithm.h>
#include <itkProgressReporter.h>
#include <limits.h>
#include "itkHelpers.h"

namespace itk {

template<class TInputImage, class TOutputImage, class TSearch,class TPatch, class TPreselectionFilter>
NonLocalPatchesFilter<TInputImage, TOutputImage, TSearch, TPatch,TPreselectionFilter>
::NonLocalPatchesFilter() : 
m_SearchKernel(),m_PatchKernel()
{
  m_DefaultBoundaryCondition.SetConstant( itk::NumericTraits<PixelType>::Zero );
  m_BoundaryCondition = &m_DefaultBoundaryCondition;
  m_PreselectionFilter = TPreselectionFilter::New();
}
  
template<class TInputImage, class TOutputImage, class TSearch,class TPatch, class TPreselectionFilter>
void 
NonLocalPatchesFilter<TInputImage, TOutputImage, TSearch, TPatch,TPreselectionFilter>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // get pointers to the input and output
  typename Superclass::InputImagePointer  inputPtr = 
    const_cast< TInputImage * >( this->GetInput() );
  
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

template<class TInputImage, class TOutputImage, class TSearch,class TPatch, class TPreselectionFilter>
void 
NonLocalPatchesFilter<TInputImage, TOutputImage, TSearch, TPatch,TPreselectionFilter>::BeforeThreadedGenerateData(void)
{
  //DEBUG
  m_Weights=TOutputImage::New();
  m_SmallestWeights=TOutputImage::New();
  
  allocate_same<TOutputImage,TOutputImage>(m_Weights,this->GetInput());
  allocate_same<TOutputImage,TOutputImage>(m_SmallestWeights,this->GetInput());
  
  m_Weights->FillBuffer(0.0);
  m_SmallestWeights->FillBuffer(0.0);
  
  //TODO: check if ROI have the same dimensions as input image (if present)

}


#if (ITK_VERSION_MAJOR==3)
template<class TInputImage, class TOutputImage, class TSearch,class TPatch, class TPreselectionFilter>
void
NonLocalPatchesFilter<TInputImage, TOutputImage, TSearch, TPatch,TPreselectionFilter>
::ThreadedGenerateData (const OutputImageRegionType& 
                              outputRegionForThread,
                              int threadId)
#else
template<class TInputImage, class TOutputImage, class TSearch,class TPatch, class TPreselectionFilter>
void
NonLocalPatchesFilter<TInputImage, TOutputImage, TSearch, TPatch,TPreselectionFilter>
::ThreadedGenerateData (const OutputImageRegionType &outputRegionForThread,
                              itk::ThreadIdType threadId)
#endif
{
  // Neighborhood iterators
  NeighborhoodIteratorType search_iter;
  NeighborhoodIteratorType patch_iter1;
  NeighborhoodIteratorType patch_iter2;

  // Find the boundary "faces"
  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList;
  
  itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> fC;
  
  faceList = fC(this->GetInput(), outputRegionForThread, m_SearchKernel.GetRadius());

  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;

  itk::ImageRegionIterator<TOutputImage> o_iter;
  itk::ImageRegionConstIterator<RoiImageType> roi_iterator;

  itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  // Process the boundary faces, these are N-d regions which border the
  // edge of the buffer

  const SearchKernelIteratorType searchKernelBegin = m_SearchKernel.Begin();
  const SearchKernelIteratorType searchKernelEnd   = m_SearchKernel.End();

  const SearchKernelIteratorType patchKernelBegin = m_PatchKernel.Begin();
  const SearchKernelIteratorType patchKernelEnd   = m_PatchKernel.End();
  //for thread safety
  PreselectionFilterPointerType preselection=m_PreselectionFilter;
  
  for (fit = faceList.begin(); fit != faceList.end(); ++fit)
  { 
    search_iter = NeighborhoodIteratorType(m_SearchKernel.GetRadius(),
                                      this->GetInput(), *fit);

    patch_iter1 = NeighborhoodIteratorType(m_PatchKernel.GetRadius(),
                                      this->GetInput(), *fit);

    patch_iter2 = NeighborhoodIteratorType(m_PatchKernel.GetRadius(),
                                      this->GetInput(), *fit);
    
    o_iter = itk::ImageRegionIterator<OutputImageType>(this->GetOutput(), *fit);
    
    search_iter.OverrideBoundaryCondition(m_BoundaryCondition);
    patch_iter1.OverrideBoundaryCondition(m_BoundaryCondition);
    patch_iter2.OverrideBoundaryCondition(m_BoundaryCondition);
    search_iter.GoToBegin();
    
    if(m_RoiImage.IsNotNull())
    {
      roi_iterator = itk::ImageRegionConstIterator<RoiImageType>(m_RoiImage,*fit);
      roi_iterator.GoToBegin();
      
      while ( ! o_iter.IsAtEnd() )
      {
        if(roi_iterator.Get()==0)
        {
          ++roi_iterator;
          ++search_iter;
          o_iter.Set(0);//TODO: maybe put input value here?
          ++o_iter;
          ++patch_iter1;
          progress.CompletedPixel();
          continue;
        } 
        ++roi_iterator;
        o_iter.Set( this->Evaluate( search_iter, patch_iter1, patch_iter2, searchKernelBegin, searchKernelEnd, patchKernelBegin, patchKernelEnd,preselection) );
        ++search_iter;
        ++o_iter;
        ++patch_iter1;
        progress.CompletedPixel();
      }
    } else {
      while ( ! o_iter.IsAtEnd() )
      {
        o_iter.Set( this->Evaluate( search_iter, patch_iter1, patch_iter2, searchKernelBegin, searchKernelEnd, patchKernelBegin, patchKernelEnd,preselection) );
        ++search_iter;
        ++o_iter;
        ++patch_iter1;
        progress.CompletedPixel();
      }
    }
  }
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch, class TPreselectionFilter>
void
NonLocalPatchesFilter<TInputImage, TOutputImage, TSearch, TPatch, TPreselectionFilter>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Search Kernel: " << m_SearchKernel << std::endl;
  os << indent << "Patch Kernel: "  << m_PatchKernel  << std::endl;
  os << indent << "Boundary condition: "  << typeid( *m_BoundaryCondition ).name() << std::endl;
  os << indent << "Preselection filter: " << typeid( m_PreselectionFilter ).name() << std::endl;
}

}// end namespace itk

#endif

// kate: space-indent on; hl c++;indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;tab-width 2 