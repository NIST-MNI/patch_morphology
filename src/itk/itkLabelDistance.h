#ifndef __mincLabelDistance_h__
#define __mincLabelDistance_h__

#include <itkOffset.h>


namespace itk
{
  //! Label distance metric
  template <class LabelPixelType,unsigned int dim>
    struct LabelDistance
  {
    typedef itk::Offset<dim> TOffset;
    
    LabelPixelType label;
    double         distance;
    size_t         sample;
    double         offset;
    double         grading;
    
    LabelDistance(LabelPixelType _label,
                  double         _distance,
                  size_t         _sample,
                  const TOffset& _offset,
                  double         _grading
                 ):
      label(_label),distance(_distance),sample(_sample),grading(_grading)
      {
        offset=0;
        for(int i=0;i<dim;i++)
          offset+=_offset[i]*_offset[i];
      }
  };
} //minc

#endif //__mincLabelDistance_h__
