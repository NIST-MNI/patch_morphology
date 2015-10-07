#ifndef __mincMedianWeightedDiscreteLabeling_h__
#define __mincMedianWeightedDiscreteLabeling_h__

#include <vector>
#include <algorithm>
#include <itkOffset.h>
#include <itkLightObject.h>

namespace itk
{
  
  template<class LabelPixelType,class OutputPixel,class TWeight,unsigned int dim>
        class MedianWeightedDiscreteLabeling: public itk::LightObject
  {
  public:
    typedef LabelDistance<LabelPixelType,dim> LabelDistanceType;
    typedef std::vector<LabelDistanceType> LabelDistanceVector;

    typedef MedianWeightedDiscreteLabeling             Self;
    itkTypeMacro(MedianWeightedDiscreteLabeling, itk::LightObject);    
    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    itkNewMacro(Self);
    
  protected:
    TWeight m_Weight;
    size_t  m_LabelCount;
    
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
    
    void SetLabelCount(size_t c)
    {
      m_LabelCount=c;
    }
    
    size_t GetLabelCount(void) const
    {
      return m_LabelCount;
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
      
      if(samples.empty())
        return; 
      
      compare_labels cmp;
      std::sort(samples.begin(),samples.end(),cmp);
      
      double median_distance=samples[samples.size()/2].distance;
      
      if(median_distance<distance_epsilon)
        median_distance=distance_epsilon;
      
      double total_weight=0;
      size_t out_label=0;
      double max_weight=0;
      
      std::vector<double> weights(m_LabelCount,0.0);
      std::vector<double> sdistances(m_LabelCount,0.0);
      std::vector<double> gradings(m_LabelCount,0.0);
      std::vector<int>    counts(m_LabelCount,0);
      
      for(::size_t i=0;i<samples.size();++i)
      {
        double w=m_Weight(samples[i].distance,median_distance);
        weights[samples[i].label]+=w;
        total_weight+=w;
        sdistances[samples[i].label]+=samples[i].distance;
        counts[samples[i].label]++;
        gradings[samples[i].label]+=samples[i].grading*w;
      }

      for(::size_t i=0;i<m_LabelCount;++i)
      {
        if(weights[i]>=max_weight) {
          max_weight=weights[i];
          out_label=i;
        }
      }
      output_value=out_label;
      confidence=max_weight/total_weight;
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
    
    bool operator!=(const MedianWeightedDiscreteLabeling &a)
    {
      return true;
    }
  };
} //minc

#endif //__mincMedianWeightedDiscreteLabeling_h__