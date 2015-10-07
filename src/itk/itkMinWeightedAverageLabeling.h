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
#ifndef __mincMinWeightedAverageLabeling_h__
#define __mincMinWeightedAverageLabeling_h__

#include <vector>
#include <algorithm>
#include <itkOffset.h>
#include <itkLightObject.h>

namespace itk
{
  template<class LabelPixelType,class OutputPixel,class TWeight,unsigned int dim>
    class MinWeightedAverageLabeling: public itk::LightObject
  {
  public:
    typedef LabelDistance<LabelPixelType,dim> LabelDistanceType;
    typedef std::vector<LabelDistanceType> LabelDistanceVector;

    typedef MinWeightedAverageLabeling             Self;
    itkTypeMacro(MinWeightedAverageLabeling, itk::LightObject);
    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    itkNewMacro(Self);
    
  protected:
    TWeight m_Weight;
  public:
    
    void SetWeight(const TWeight& w)
    {
      m_Weight=w;
    }
    
    TWeight& GetWeight(void)
    {
      return m_Weight;
    }
    
    OutputPixel default_value(void) const
    {
      return 0;
    }
    
    size_t GetLabelCount(void) const
    {
      return 2;
    }
    
    void process(LabelDistanceVector &samples,
                    OutputPixel &output_value,
                    double &confidence,
                    double &search_distance,
                    double &grading,
                    std::vector<double>& prob
                ) const
    {
      output_value=default_value();
      confidence=0;
      search_distance=0;
      grading=0;
      
      const double distance_epsilon=1e-4; //compatible with Pierrick's code
      
      //TODO: make up a meaningful value
      if(samples.empty())
        return;
      
      double min_distance=samples[0].distance;
      
      // search for the minimal value
      for(::size_t i=1;i<samples.size();++i)
        if(min_distance>samples[i].distance)
          min_distance=samples[i].distance;
      
      //make sure that minimal distance is not 0
      if(min_distance<distance_epsilon)
        min_distance=distance_epsilon;
      
      double out=0;
      double total_weight=0;
      
      for(::size_t i=0;i<samples.size();++i)
      {
        double w=m_Weight(samples[i].distance,min_distance);
        double v=w*samples[i].label;
        out+=v;
        total_weight+=w;
        //search_distance+=v*sqrt(samples[i].offset);
        search_distance+=samples[i].distance;
        grading+=w*samples[i].grading;
      }
      
      output_value=out/total_weight;
      confidence=std::max((double)output_value,1.0-output_value);  //it is either one ore another
      
      //search_distance/=total_weight;
      search_distance/=samples.size(); 
      
      grading/=total_weight;

      if(prob.size()==2)
      {
        prob[0]=1.0-output_value;
        prob[1]=output_value;
      }
    }
    
    bool operator!=(const MinWeightedAverageLabeling &a)
    {
      return true;
    }
  };

  template<class LabelPixelType,class OutputPixel,class TWeight,unsigned int dim>
      std::ostream & operator<<(std::ostream &os, const MinWeightedAverageLabeling<LabelPixelType,OutputPixel,TWeight,dim> &w)
  {
    os << "[ MinWeightedAverageLabeling<LabelPixelType,OutputPixel,TWeight> ]";
    return os;
  }

  
} //minc

#endif //__mincMinWeightedAverageLabeling_h__