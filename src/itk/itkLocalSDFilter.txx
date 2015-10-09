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
#ifndef _MINCLOCALSDFILTER_TXX
#define _MINCLOCALSDFILTER_TXX

namespace itk {

template<class TInputImage, class TOutputImage, class TKernel>
LocalSDFilter<TInputImage, TOutputImage, TKernel>
::LocalSDFilter()
{
}

template<class TInputImage, class TOutputImage, class TKernel>
typename LocalSDFilter<TInputImage, TOutputImage, TKernel>::PixelType
LocalSDFilter<TInputImage, TOutputImage, TKernel>
::Evaluate(const NeighborhoodIteratorType &nit,
          const KernelIteratorType kernelBegin,
          const KernelIteratorType kernelEnd)
{
  unsigned int i;
  double sum=0.0,sum2=0.0;
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
      double v=nit.GetPixel(i);
      sum  += v;
      sum2 += v*v;
      cnt++;
      }
    }
  sum/=cnt;
  sum2/=cnt;
  return sum2-sum*sum;
}


}// end namespace itk


#endif //_MINCLOCALSDFILTER_TXX

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2;hl c++

