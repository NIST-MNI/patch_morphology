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
#ifndef _mincVectorL2PatchDistance_h_
#define _mincVectorL2PatchDistance_h_

#include <iostream>
#include <vector>
#include <algorithm>
#include <itkLightObject.h>

namespace itk
{ 
  template<class T1,class T2> T1 WeightedL2Distance(
      const typename itk::VariableLengthVector<T2>& a,
      const typename itk::VariableLengthVector<T2>& b,
      const typename itk::VariableLengthVector<T2>& w)
  {
  
    T1 _sum = 0.0;

    for ( unsigned int i = 0; i < a.GetSize(); i++ )
    {
      T1 _d= static_cast<T1>(a[i]) - static_cast<T1>(b[i]);
      T1 _w= static_cast<T1>(w[i]);
      _sum += _d*_d*_w;
    }
  
    return _sum;
  }
  
  template<class TImage,class TPatch> 
    class VectorL2PatchDistance : public itk::LightObject
  {
  public:
    typedef VectorL2PatchDistance Self;
    itkTypeMacro(VectorL2PatchDistance, itk::LightObject);
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    itkNewMacro(Self);
    typedef typename TImage::InternalPixelType PixelType;
    typedef typename itk::VariableLengthVector<PixelType> VectorType;


  protected:
  /** Neighborhood iterator type. */
    typedef typename itk::ConstNeighborhoodIterator<TImage>  NeighborhoodIteratorType;
    typedef typename TPatch::ConstIterator  PatchKernelIteratorType;
    typedef typename TPatch::PixelType      KernelPixelType;
    
    VectorType m_Weights;
  public:
    double distance(const NeighborhoodIteratorType &patchIt1,
                    const NeighborhoodIteratorType &patchIt2,
                    const PatchKernelIteratorType patchKernelBegin,
                    const PatchKernelIteratorType patchKernelEnd) const
    {
      ::size_t i;
      PatchKernelIteratorType kernel_it;
      ::size_t center = (::size_t) (patchIt1.Size() / 2); // get offset of center pixel

      double total=0.0;
      int cnt=0;
      
      for( i=0, kernel_it=patchKernelBegin; kernel_it<patchKernelEnd; ++kernel_it, ++i )
      {
        // if structuring element is positive, use the pixel under that element
        // in the image

        if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero  )
        {
          total+=
            WeightedL2Distance<double,PixelType>(
              patchIt1.GetPixel( i ),
              patchIt2.GetPixel( i ),
              m_Weights );
          cnt++;
        }
      }
      return total/cnt;
    }

    bool operator!=(const Self &a)
    {
      return false;
    }
    
    void set_weights(const std::vector<double> & _w)
    {
      m_Weights.SetSize(_w.size());
      for(size_t i=0;i<_w.size();i++)
        m_Weights[i]=_w[i];
    }
  };
}//minc

#endif //_mincVectorL2PatchDistance_h_