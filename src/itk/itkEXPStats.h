#ifndef _ITK_EXP_STATS_h_
#define _ITK_EXP_STATS_h_

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
  
  template<class TImage,class TPatch,class TDistance=L2Distance> class ExpSimilarityStats: public itk::LightObject
  {
  public:
    typedef ExpSimilarityStats Self;
    itkTypeMacro(ExpSimilarityStats, itk::LightObject);
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    itkNewMacro(Self);

  protected:
  /** Neighborhood iterator type. */
    typedef typename itk::ConstNeighborhoodIterator<TImage>  NeighborhoodIteratorType;
    typedef typename TPatch::ConstIterator  PatchKernelIteratorType;
    typedef typename TPatch::PixelType      KernelPixelType;
    ExpSimilarityStats(const Self&); //purposely not implemented

    std::vector<double> _A;
    std::vector<double> _b;
    std::vector<double> _x;
    std::vector<int>    _sample_count;
    std::vector<int>    _label_hist;
    std::vector<double> _threshold;
    
    size_t              _library_size;
    size_t              _sample_size;
    size_t              _feature_size;
    size_t              _dist_cnt;
    
    int                 _most_common_label;
    int                 _second_most_common_label;
    TDistance           _distance_op;
    int                 _cur_label;
    double              _dist_threshold;
    double              _inner_dist_threshold;

    ExpSimilarityStats():
      _cur_label(-1)
    {
    }

    ~ExpSimilarityStats()
    {
      cleanup();
    }

    void cleanup(void)
    {
    }

  public:
    Pointer clone(void) const
    {
      Pointer _clone=Self::New();
      
      _clone->_distance_op=_distance_op;
      return _clone;
    }

    void allocate_memory(size_t library_size,
                         size_t sample_size,
                         size_t feature_size,
                         size_t dist_cnt,
                         double distance_threshold,
                         double inner_distance_threshold )
    {
      cleanup();
      //need to copy for a potential multi-threading operation

      _sample_size=sample_size;
      _library_size=library_size;
      _feature_size=feature_size;
      _dist_cnt=dist_cnt;

      _A.resize(_sample_size*_library_size*feature_size);
      _b.resize(feature_size);
      _x.resize(_sample_size*_library_size);
      _sample_count.resize(_library_size);
      _threshold.resize(_library_size);

      _dist_threshold=distance_threshold;
      _inner_dist_threshold=inner_distance_threshold;
    }

    void start_new_regression(void)
    {
      //_sample=0;
      _cur_label=0;
      _sample_count.assign(_library_size,0);
      _threshold.assign(_library_size,0.0);
    }

    void set_input(const NeighborhoodIteratorType &patchIt1,
                   int _label,
                   const PatchKernelIteratorType patchKernelBegin,
                   const PatchKernelIteratorType patchKernelEnd
                  )
    {
      //double norm_k=1.0/(norm_max-norm_min);
      ::size_t i,j=0;
      PatchKernelIteratorType kernel_it;
      //::size_t center = (::size_t) (patchIt1.Size() / 2); // get offset of center pixel
      _cur_label=_label;
      for( i=0, kernel_it=patchKernelBegin; kernel_it<patchKernelEnd; ++kernel_it, ++i )
      {
        // if structuring element is positive, use the pixel under that element
        // in the image
        if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero  )
        {
          _b[j++] = patchIt1.GetPixel(i);
        }
      }
    }
    
    void add_training_sample(
                            int    label,
                            double grading,
                            const  NeighborhoodIteratorType &patchIt2,
                            const  PatchKernelIteratorType patchKernelBegin,
                            const  PatchKernelIteratorType patchKernelEnd,
                            size_t sample_id 
                            )
    {
      //TODO: add another training sample
      ::size_t i=0;
      PatchKernelIteratorType kernel_it;
      
      /*If current label is set and doesn't matchup - filter*/
      if(_cur_label>=0 && 
         label != _cur_label)
        return;
      
      //TODO: maybe calculate distances here?
      double *__A= &_A[ (_sample_count[sample_id]+sample_id*_sample_size)*_feature_size ];
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
      _sample_count[sample_id]++;
    }

    bool regress(double &threshold , double &confidence, std::vector<double> &prob, double &grading)
    {
      size_t usefull_samples=0;
      for( size_t k=0;k<_library_size;k++ )
      {
        for(size_t t=0;t<_sample_count[k];t++)
        {
          double dist=_distance_op(&_b[0],&_A[(k*_sample_size+t)*_feature_size],_feature_size);
          _x[k*_sample_size+t]=dist;
        }
        
        if( _sample_count[k]>0 )
        {
          std::sort(_x.begin() + k*_sample_size , _x.begin() + k*_sample_size + _sample_count[k]);
          _threshold[usefull_samples]=_x[k*_sample_size + _sample_count[k]*_inner_dist_threshold ]; 
          usefull_samples++;
        }
      }
      std::sort(_threshold.begin(),_threshold.begin()+usefull_samples);
      threshold=_threshold[ usefull_samples*_dist_threshold ];
      confidence=usefull_samples;
      
      //TODO:calculate more detailed stats, put into prob
      return true;
    }

    bool operator!=(const Self &a)
    {
      return false;
    }

    double default_value(void) const
    {
      return -1.0;
    }
  };

}//minc

#endif //_ITK_EXP_STATS_h_