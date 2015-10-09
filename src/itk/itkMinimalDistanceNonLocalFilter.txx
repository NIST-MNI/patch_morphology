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
#ifndef __mincMinimalDistanceNonLocalFilter_txx
#define __mincMinimalDistanceNonLocalFilter_txx


namespace itk {

template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance,class TPreselectionFilter>
MinimalDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance,TPreselectionFilter>
::MinimalDistanceNonLocalFilter()
{
  m_Distance = TDistance::New();
}

template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance,class TPreselectionFilter>
typename MinimalDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance,TPreselectionFilter>::PixelType
MinimalDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance,TPreselectionFilter>
::Evaluate(const NeighborhoodIteratorType &searchIt,
          const NeighborhoodIteratorType &patchIt1,
          NeighborhoodIteratorType &patchIt2,
          const SearchKernelIteratorType searchKernelBegin,
          const SearchKernelIteratorType searchKernelEnd,
          const PatchKernelIteratorType patchKernelBegin,
          const PatchKernelIteratorType patchKernelEnd,
          PreselectionFilterPointerType flt
          )
{
  ::size_t i;
  SearchKernelIteratorType kernel_it;
  ::size_t center = (::size_t) (searchIt.Size() / 2); // get offset of center pixel

  double smallest_distance=-1.0;
  
  for( i=0, kernel_it=searchKernelBegin; kernel_it<searchKernelEnd; ++kernel_it, ++i )
  {
    // if structuring element is positive, use the pixel under that element
    // in the image
    if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero )
    {
      patchIt2.SetLocation(searchIt.GetIndex(i)); //move patch
      
      if(!patchIt2.InBounds()) continue;
      //check preselection
      if(!flt->select(patchIt1.GetIndex(),patchIt2.GetIndex())) continue;
      
      if(i==center) 
        continue;
      
      double distance=m_Distance->distance(patchIt1,patchIt2,patchKernelBegin,patchKernelEnd);
      
      if(distance<smallest_distance || smallest_distance<0.0) 
        smallest_distance=distance;
    }
  }
  
  return (PixelType)(smallest_distance>0.0?smallest_distance:0);
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TPreselectionFilter>
void MinimalDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance,TPreselectionFilter>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


}// end namespace itk
#endif

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;show-tabs on;tab-width 2;hl c++
