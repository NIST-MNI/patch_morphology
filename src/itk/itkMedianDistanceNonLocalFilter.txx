#ifndef __mincMedianDistanceNonLocalFilter_txx
#define __mincMedianDistanceNonLocalFilter_txx

#include <vector>
#include <algorithm>

namespace itk {

template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance,class TPreselectionFilter>
MedianDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance,TPreselectionFilter>
::MedianDistanceNonLocalFilter()
{
  m_Distance = TDistance::New();
}

template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance,class TPreselectionFilter>
typename MedianDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance,TPreselectionFilter>::PixelType
MedianDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance,TPreselectionFilter>
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
  std::vector<double> dist;
  dist.reserve(searchIt.Size());

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
      
      dist.push_back(m_Distance->distance(patchIt1,patchIt2,patchKernelBegin,patchKernelEnd));
      
    }
  }
  std::sort(dist.begin(),dist.end());

  return (PixelType)(dist.size()%2?(dist[dist.size()/2]+dist[dist.size()/2+1])/2:dist[dist.size()/2]);
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TPreselectionFilter>
void MedianDistanceNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance,TPreselectionFilter>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


}// end namespace itk
#endif

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;show-tabs on;tab-width 2;hl c++
