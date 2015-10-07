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
#ifndef __mincMinWeightedDiscreteLabeling_h__
#define __mincMinWeightedDiscreteLabeling_h__

#include <vector>
#include <algorithm>
#include <itkOffset.h>
#include <itkLightObject.h>

namespace itk
{
  template<class LabelPixelType,class OutputPixel,class TWeight,unsigned int dim>
    class MinWeightedDiscreteLabeling: public itk::LightObject
  {
  public:
    typedef LabelDistance<LabelPixelType,dim>  LabelDistanceType;
    typedef std::vector<LabelDistanceType> LabelDistanceVector;

    typedef MinWeightedDiscreteLabeling             Self;
    itkTypeMacro(MinWeightedDiscreteLabeling, itk::LightObject);    
    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    itkNewMacro(Self);
    
  protected:
    TWeight m_Weight;
    size_t  m_LabelCount;
  public:
    
    void SetWeight(const TWeight& w)
    {
      m_Weight=w;
    }
    
    TWeight& GetWeight(void)
    {
      return m_Weight;
    }
    
    void SetLabelCount(size_t c)
    {
      m_LabelCount=c;
    }
    
    size_t GetLabelCount(void) const
    {
      return m_LabelCount;
    }
    
    OutputPixel default_value(void) const
    {
      return 0;
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
        return ; 
      
      double min_distance=samples[0].distance;
      
      for(::size_t i=1;i<samples.size();++i)
        if(min_distance>samples[i].distance) min_distance=samples[i].distance;

      if(min_distance<distance_epsilon)
        min_distance=distance_epsilon;

      double total_weight=0;
      double max_weight=0.0;
      size_t out_label=0;
      
      std::vector<double> weights(m_LabelCount,0.0);
      std::vector<double> sdistances(m_LabelCount,0.0);
      std::vector<double> gradings(m_LabelCount,0.0);
      std::vector<int>    counts(m_LabelCount,0);
      
      for(::size_t i=0;i<samples.size();++i)
      {
        if(samples[i].label >= m_LabelCount) abort();
        
        double w=m_Weight(samples[i].distance,min_distance);
        weights[samples[i].label]+=w;
        sdistances[samples[i].label]+=samples[i].distance;
        counts[samples[i].label]++;
        total_weight+=w;
        gradings[samples[i].label]+=w*samples[i].grading;
      }
      
      for(::size_t i=0;i<m_LabelCount;++i)
      {
        if(weights[i]>=max_weight) {
          max_weight=weights[i];
          out_label=i;
          //out_label=weights[i]*100/total_weight;
        }
      }
      confidence=max_weight/total_weight; 
      output_value=out_label;
      if(counts[out_label]>0)
        search_distance=sdistances[out_label]/counts[out_label];
      else
        search_distance=0.0;
      
      grading=gradings[out_label]/total_weight;
      
      if(prob.size()==m_LabelCount)
      {
        for(::size_t i=0;i<m_LabelCount;++i)
        {
          prob[i]=weights[i]/total_weight;
        }
      }
    }
    
    bool operator!=(const MinWeightedDiscreteLabeling &a)
    {
      return true;
    }
  };

} //minc

#endif //__mincMinWeightedDiscreteLabeling_h__