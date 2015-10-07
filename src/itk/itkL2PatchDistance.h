#ifndef _mincL2PatchDistance_h_
#define _mincL2PatchDistance_h_

#include <iostream>
#include <vector>
#include <algorithm>
#include <itkLightObject.h>

namespace itk
{
  template<class TImage,class TPatch> class L2PatchDistance : public itk::LightObject
  {
  public:
    typedef L2PatchDistance Self;
    itkTypeMacro(L2PatchDistance, itk::LightObject);    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    itkNewMacro(Self);


  protected:
  /** Neighborhood iterator type. */
    typedef typename itk::ConstNeighborhoodIterator<TImage>  NeighborhoodIteratorType;
    typedef typename TPatch::ConstIterator  PatchKernelIteratorType;
    typedef typename TPatch::PixelType      KernelPixelType;

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
          double v=patchIt1.GetPixel(i)-patchIt2.GetPixel(i);
          total+=v*v;
          cnt++;
        }
      }
      return total/cnt;
    }

    bool operator!=(const Self &a)
    {
      return false;
    }
  };

}//minc

#endif //_mincL2PatchDistance_h_