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
#ifndef __mincLocalMeanFilter_txx__
#define __mincLocalMeanFilter_txx__

namespace itk {

template<class TInputImage, class TOutputImage, class TKernel>
LocalMeanFilter<TInputImage, TOutputImage, TKernel>
::LocalMeanFilter()
{
}

template<class TInputImage, class TOutputImage, class TKernel>
typename LocalMeanFilter<TInputImage, TOutputImage, TKernel>::PixelType
LocalMeanFilter<TInputImage, TOutputImage, TKernel>
::Evaluate(const NeighborhoodIteratorType &nit,
          const KernelIteratorType kernelBegin,
          const KernelIteratorType kernelEnd)
{
  unsigned int i;
  double sum=0.0;
  unsigned int cnt=0;

  KernelIteratorType kernel_it;

  for( i=0, kernel_it=kernelBegin; kernel_it<kernelEnd; ++kernel_it, ++i )
    {
    // if structuring element is positive, use the pixel under that element
    // in the image
    if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero )
      {
      // note we use GetPixel() on the SmartNeighborhoodIterator to
      // respect boundary conditions
      sum += nit.GetPixel(i);
      cnt++;
      }
    }
  
  return sum/cnt;
}


}// end namespace itk


#endif //__mincLocalMeanFilter_txx__

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl C++
