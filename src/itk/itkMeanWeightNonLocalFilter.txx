#ifndef __mincMeanWeightNonLocalFilter_txx
#define __mincMeanWeightNonLocalFilter_txx


namespace itk {

template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
MeanWeightNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::MeanWeightNonLocalFilter():
  m_sigma2(1.0)
{
  m_Distance = TDistance::New();
}

template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
typename MeanWeightNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>::PixelType
MeanWeightNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
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

	double total_weight=0.0;
	double smallest_weight=0.0;
	
  for( i=0, kernel_it=searchKernelBegin; kernel_it<searchKernelEnd; ++kernel_it, ++i )
	{
    // if structuring element is positive, use the pixel under that element
    // in the image
    if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero )
		{
			patchIt2.SetLocation(searchIt.GetIndex(i)); //move patch
			
			if(!patchIt2.InBounds()) continue;
			
			if(!flt->select(patchIt1.GetIndex(),patchIt2.GetIndex())) continue;

			//TODO: move distance & weight calculations into templates
			double distance=0;
			
			if(i!=center) 
				distance=m_Distance->distance(patchIt1,patchIt2,patchKernelBegin,patchKernelEnd);

			double weight=m_Weight(distance,m_sigma2);
      total_weight+=weight;
			if(i!=center && (weight>smallest_weight || smallest_weight==0.0))
				smallest_weight=weight;
		}
	}
	//add central voxel
	//DEBUG
	this->m_Weights->SetPixel(searchIt.GetIndex(center),total_weight);
	this->m_SmallestWeights->SetPixel(searchIt.GetIndex(center),smallest_weight);
	
  //return average weight
  return (PixelType)(total_weight/i);
}

template<class TInputImage, class TOutputImage, class TSearch,class TPatch,class TDistance,class TWeight,class TPreselectionFilter>
void MeanWeightNonLocalFilter<TInputImage, TOutputImage, TSearch, TPatch, TDistance, TWeight,TPreselectionFilter>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "sigma^2: " << m_sigma2 << std::endl;
}

}// end namespace itk
#endif

// kate: space-indent on; hl c++;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 2 