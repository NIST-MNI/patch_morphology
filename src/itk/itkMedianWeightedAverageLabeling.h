#ifndef __mincMedianWeightedAverageLabeling_h__
#define __mincMedianWeightedAverageLabeling_h__

#include <vector>
#include <algorithm>
#include <itkOffset.h>
#include <itkLightObject.h>

namespace itk
{
  template<class LabelPixelType,class OutputPixel,class TWeight,unsigned int dim>
    class MedianWeightedAverageLabeling: public itk::LightObject
  {
  public:
    typedef LabelDistance<LabelPixelType,dim> LabelDistanceType;
    typedef std::vector<LabelDistanceType> LabelDistanceVector;

    typedef MedianWeightedAverageLabeling             Self;
    itkTypeMacro(MedianWeightedAverageLabeling, itk::LightObject);
    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    itkNewMacro(Self);
    
  protected:
    TWeight m_Weight;
    
    struct compare_labels
    {
      bool operator()(const LabelDistanceType &a,const LabelDistanceType &b) const
      {
        return a.distance<b.distance;
      }
    };
    
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
    
    void process(LabelDistanceVector &samples,
                    OutputPixel &output_value,
                    double &confidence,
                    double &search_distance,
                    double &grading,
                    std::vector<double>& prob
                ) const
    {
      output_value=default_value();
      confidence=0.0;
      search_distance=0.0;
      grading=0.0;
      
      const double distance_epsilon=1e-4; //compatible with Pierrick's code
      
      if(samples.empty())
        return ; 
      
      compare_labels cmp;
      std::sort(samples.begin(),samples.end(),cmp);
      
      double median_distance=samples[samples.size()/2].distance;
      
        
      if(median_distance<distance_epsilon)
        median_distance=distance_epsilon;
      
      double out=0;
      double total_weight=0;
      
      for(::size_t i=0;i<samples.size();++i)
      {
        double w=m_Weight(samples[i].distance,median_distance);
        double v=w*samples[i].label;
        out+=v;
        total_weight+=w;
        search_distance+=samples[i].distance;
        grading+=w*samples[i].grading;
      }
      
      output_value=out/total_weight;
      confidence=std::max((double)output_value,1.0-output_value);  //it is either one ore another
      search_distance/=samples.size();
      grading/=total_weight;
      
      if(prob.size()==2)
      {
        prob[0]=1.0-output_value;
        prob[1]=output_value;
      }
    }
    
    bool operator!=(const MedianWeightedAverageLabeling &a)
    {
      return true;
    }
  };

} //minc

#endif //__mincMedianWeightedAverageLabeling_h__
