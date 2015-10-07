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
