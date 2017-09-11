#ifndef _ITK_NNLS_REGRESSION_GRADING_h_
#define _ITK_NNLS_REGRESSION_GRADING_h_

#include <iostream>
#include <vector>
#include <algorithm>
#include <itkLightObject.h>

#ifdef HAVE_OPENBLAS
#include "nnls.h"
#else
#define NNLS_HANDLE void *
#endif

namespace itk
{
  //WIP: grading is not working yet!
  template<class TImage,class TPatch> class NNLSRegressionGrading: public itk::LightObject
  {
  public:
    typedef NNLSRegressionGrading Self;
    itkTypeMacro(L1PatchDistance, itk::LightObject);    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    itkNewMacro(Self);
    
    double        norm_min,norm_max;
    
    
  protected:
  /** Neighborhood iterator type. */
    typedef typename itk::ConstNeighborhoodIterator<TImage>  NeighborhoodIteratorType;
    typedef typename TPatch::ConstIterator  PatchKernelIteratorType;
    typedef typename TPatch::PixelType      KernelPixelType;
    NNLSRegressionGrading(const Self&); //purposely not implemented
    
    std::vector<double> NNLS_A;
    std::vector<double> NNLS_b;
    std::vector<double> NNLS_x;
    std::vector<double> label_weight;
    std::vector<int>    labels;
    int                 nnls_sample;
    size_t              _sample_size,_feature_size,_label_cnt;
    NNLS_HANDLE         NNLS_work;
    
    
    NNLSRegressionGrading():
      norm_min(0.0),
      norm_max(1.0),
      NNLS_work(NULL)
    {
    }
    
    ~NNLSRegressionGrading()
    {
      cleanup();
    }
    
    void cleanup(void)
    {
#ifdef HAVE_OPENBLAS      
      if(NNLS_work)
        free_nnls(NNLS_work);
#endif      
      NNLS_work=NULL;
    }

  public:
    Pointer clone(void) const
    {
      Pointer _clone=Self::New();
      _clone->norm_min=norm_min;
      _clone->norm_max=norm_max;
      
      return _clone;
    }
    
    void allocate_memory(size_t sample_size, size_t feature_size,size_t label_cnt)
    {
      cleanup();
      _sample_size=sample_size;
      _feature_size=feature_size;
      _label_cnt=label_cnt;
//       std::cout<<"_feature_size="<<_feature_size<<"\t"
//                <<"_sample_size="<<_sample_size<<"\t"
//                <<"_label_cnt="<<_label_cnt<<std::endl;
               
      NNLS_A.resize(sample_size*feature_size);
      NNLS_b.resize(feature_size);
      NNLS_x.resize(sample_size);
      
      labels.resize(sample_size);
      
      label_weight.resize(label_cnt);
#ifdef HAVE_OPENBLAS      
      NNLS_work=allocate_nnls(1,_feature_size, _sample_size);
#else
      fprintf(stderr,"Compiled without OpenBLAS, NNLS not implemented\n");
      abort();
#endif
    }
    
    
    void  start_new_regression(void)
    {
      nnls_sample=0;
    }
    
    
    void set_input(const NeighborhoodIteratorType &patchIt1,
                   const PatchKernelIteratorType patchKernelBegin,
                   const PatchKernelIteratorType patchKernelEnd
                  )
    {
      double norm_k=1.0/(norm_max-norm_min);
      ::size_t i,j=0;
      PatchKernelIteratorType kernel_it;
      //::size_t center = (::size_t) (patchIt1.Size() / 2); // get offset of center pixel
      
      for( i=0, kernel_it=patchKernelBegin; kernel_it<patchKernelEnd; ++kernel_it, ++i )
      {
        // if structuring element is positive, use the pixel under that element
        // in the image
        if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero  )
        {
          NNLS_b[j++]=(patchIt1.GetPixel(i)-norm_min)*norm_k;
        }
      }
    }
    
    void add_training_sample(
                            int label,
                            double grading,
                            const NeighborhoodIteratorType &patchIt2,
                            const PatchKernelIteratorType patchKernelBegin,
                            const PatchKernelIteratorType patchKernelEnd,
                            size_t sample_id 
                            )
    {
      //TODO: add another training sample
      double norm_k=1.0/(norm_max-norm_min);
      ::size_t i=0,j=0;
      PatchKernelIteratorType kernel_it;
      
      if(nnls_sample>=_sample_size) abort();
      
      labels[nnls_sample]=label;
      
      for( kernel_it=patchKernelBegin; kernel_it<patchKernelEnd; ++kernel_it, ++i )
      {
        // if structuring element is positive, use the pixel under that element
        // in the image
        
        if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero  )
        {
          NNLS_A[ j+ nnls_sample*_feature_size ]=(patchIt2.GetPixel(i)-norm_min)*norm_k;
          j++;
        }
      }
      nnls_sample++;
    }

    
    bool regress(int &label,double &confidence,std::vector<double> &prob,double &grading)
    {
#ifdef HAVE_OPENBLAS
      nnls2_updates_single(NNLS_work,&NNLS_A[0], &NNLS_b[0], &NNLS_x[0], 0, 
                   (_feature_size+_sample_size)*2, (_feature_size+_sample_size)*2, 1,  
                   _feature_size, _sample_size, 1e-20);      
#endif

      label_weight.assign(_label_cnt,0.0);
      double total_weight=0.0;

      for(size_t k=0;k<_sample_size;k++)
      {
        total_weight+=NNLS_x[k];
        label_weight[labels[k]]+=NNLS_x[k];
      }
      double max_weight=0.0;
      size_t chosen_label=-1;
      
      for(size_t k=0;k<_label_cnt;k++)
      {
        label_weight[k]/=total_weight;
        if(label_weight[k]>max_weight)
        {
          max_weight=label_weight[k];
          chosen_label=k;
        }
        prob[k]=label_weight[k];
      }
      
      label=chosen_label;
      grading = total_weight;
      
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
  

}//itk

#endif //_ITK_NNLS_REGRESSION_GRADING_h_