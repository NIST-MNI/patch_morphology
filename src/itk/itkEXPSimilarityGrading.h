#ifndef _ITK_EXP_SIMILARITY_LABELLING_h_
#define _ITK_EXP_SIMILARITY_LABELLING_h_

#include <iostream>
#include <vector>
#include <algorithm>
#include <itkLightObject.h>

namespace itk
{

  class L2Distance
  {
  public:
    double operator() (double *a,double *b,int cnt) const
    {
      double sum=0.0;
      for(int i=0;i<cnt;i++)
        sum+=(a[i]-b[i])*(a[i]-b[i]);
      return sum/cnt;
    }
  };
  
  class L1Distance
  {
  public:
    double operator() (double *a,double *b,int cnt) const
    {
      double sum=0.0;
      for(int i=0;i<cnt;i++)
        sum+=fabs(a[i]-b[i]);
      return sum/cnt;
    }
  };
  
  template<class TImage,class TPatch,class TDistance=L2Distance> class ExpSimilarityGrading: public itk::LightObject
  {
  public:
    typedef ExpSimilarityGrading Self;
    itkTypeMacro(L1PatchDistance, itk::LightObject);
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    itkNewMacro(Self);

  protected:
  /** Neighborhood iterator type. */
    typedef typename itk::ConstNeighborhoodIterator<TImage>  NeighborhoodIteratorType;
    typedef typename TPatch::ConstIterator  PatchKernelIteratorType;
    typedef typename TPatch::PixelType      KernelPixelType;
    ExpSimilarityGrading(const Self&); //purposely not implemented
    
    std::vector<double> _A;
    std::vector<double> _b;
    std::vector<double> _x;
    std::vector<double> _label_weight;
    std::vector<double> _label_grade;
    std::vector<int>    _labels;
    std::vector<int>    _label_hist;
    std::vector<double> _grades;
    
    int                 _sample;
    size_t              _sample_size,_feature_size,_label_cnt;
    int                 _most_common_label;
    int                 _second_most_common_label;
    TDistance           _distance_op;
    double              _beta;
    
    ExpSimilarityGrading():
        _beta(0.5)
    {
    }
    
    ~ExpSimilarityGrading()
    {
      cleanup();
    }
    
    void cleanup(void)
    {
    }

  public:
      
    void set_beta(double b)
    {
        _beta=b;
    }
    
    double get_beta(void) const
    {
        return _beta;
    }
      
    Pointer clone(void) const
    {
      Pointer _clone=Self::New();
      _clone->_beta=this->_beta;
      
      return _clone;
    }
    
    void allocate_memory(size_t sample_size, 
                         size_t feature_size,
                         size_t label_cnt, 
                         TDistance distance_op=TDistance())
    {
      cleanup();
      //need to copy for a potential multi-threading operation
      _distance_op=distance_op;
      
      _sample_size=sample_size;
      _feature_size=feature_size;
      _label_cnt=label_cnt;
               
      _A.resize(sample_size*feature_size);
      _b.resize(feature_size);
      _x.resize(sample_size);
      
      _labels.resize(sample_size);
      _grades.resize(sample_size);
      
      _label_weight.resize(label_cnt);
      _label_grade.resize(label_cnt);
      _label_hist.resize(label_cnt);
    }
    
    void start_new_regression(void)
    {
      _sample=0;
      _label_hist.assign(_label_cnt,0);
      _most_common_label=-1;
      _second_most_common_label=-1;
    }
    
    
    void set_input(const NeighborhoodIteratorType &patchIt1,
                   const PatchKernelIteratorType patchKernelBegin,
                   const PatchKernelIteratorType patchKernelEnd
                  )
    {
      //double norm_k=1.0/(norm_max-norm_min);
      ::size_t i,j=0;
      PatchKernelIteratorType kernel_it;
      //::size_t center = (::size_t) (patchIt1.Size() / 2); // get offset of center pixel
      
      for( i=0, kernel_it=patchKernelBegin; kernel_it<patchKernelEnd; ++kernel_it, ++i )
      {
        // if structuring element is positive, use the pixel under that element
        // in the image
        if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero  )
        {
          _b[j++]=patchIt1.GetPixel(i);
        }
      }
    }
    
    void add_training_sample(
                            int    label,
                            double grading,
                            const NeighborhoodIteratorType &patchIt2,
                            const PatchKernelIteratorType patchKernelBegin,
                            const PatchKernelIteratorType patchKernelEnd,
                            size_t sample_id 
                            )
    {

      
      //TODO: add another training sample
      //double norm_k=1.0/(norm_max-norm_min);
      ::size_t i=0;
      PatchKernelIteratorType kernel_it;
      
      //THIS things should not happen
//       if(_sample>=_sample_size) abort();
//       
//       if(label<0          ) abort();
//       if(label>=_label_cnt) abort();
//       _label_hist[label]++;
      
      _labels[_sample]=label;
      _grades[_sample]=grading;
      
      double *__A=&_A[ _sample*_feature_size ];
      for( kernel_it=patchKernelBegin; kernel_it<patchKernelEnd; ++kernel_it, ++i )
      {
        // if structuring element is positive, use the pixel under that element
        // in the image
        if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero  )
        {
          *__A=patchIt2.GetPixel(i);
          __A++;
        }
      }
      _sample++;
    }

    bool regress(int &label , double &confidence, std::vector<double> &prob, double &grading)
    {
      const double _epsilon=1e-4;      
      
      _label_weight.assign(_label_cnt,0.0);
      _label_grade.assign(_label_cnt,0.0);
      double total_weight=0.0;
      double min_distance=1e10;
      //double avg_dist=0.0;
      
      for(size_t k=0;k<_sample_size;k++)
      {
        double dist=_distance_op(&_b[0],&_A[k*_feature_size],_feature_size);
        //avg_dist+=dist;
        _x[k]=dist;
        if(min_distance>dist || k==0) min_distance=dist;
      }
      //avg_dist/=_sample_size;
      
      if(min_distance<_epsilon) min_distance=_epsilon;
      
      //TODO: make inverse exponent another template class?
      for(size_t k=0;k<_sample_size;k++)
      {
        double w=exp(-_x[k]/(2.0*min_distance*_beta));
        total_weight+=w;
        _label_weight[_labels[k]]+=w;
        _label_grade[_labels[k]]+=w*_grades[k];
      }
      
      double max_weight=0.0;
      int chosen_label=-1;

      for(size_t k=0;k<_label_cnt;k++)
      {
        _label_weight[k]/=total_weight;
        _label_grade[k]/=total_weight;
        
        if(_label_weight[k]>max_weight || k==0)
        {
          max_weight=_label_weight[k];
          chosen_label=k;
        }
        prob[k]=_label_weight[k];
      }
      //assign grading of the most probable label 
      
      grading =_label_grade[chosen_label];
      label   = chosen_label;
      return true;
    }

    bool operator!=(const Self &a)
    {
      return false;
    }

    int default_value(void) const
    {
      return 0;
    }
  };

}//minc

#endif //_ITK_EXP_SIMILARITY_LABELLING_h_