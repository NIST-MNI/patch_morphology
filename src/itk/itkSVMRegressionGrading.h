#ifndef _ITK_SVM_REGRESSION_GRADING_h_
#define _ITK_SVM_REGRESSION_GRADING_h_

#include <iostream>
#include <vector>
#include <algorithm>
#include <itkLightObject.h>
#include "svm.h"

namespace itk
{
  //WIP: grading is not working yet!
  template<class TImage,class TPatch> class SVMRegressionGrading: public itk::LightObject
  {
  public:
    typedef SVMRegressionGrading Self;
    itkTypeMacro(L1PatchDistance, itk::LightObject);    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    itkNewMacro(Self);
    
    svm_parameter svm_param;
    double        norm_min,norm_max;
    
  protected:
  /** Neighborhood iterator type. */
    typedef typename itk::ConstNeighborhoodIterator<TImage>  NeighborhoodIteratorType;
    typedef typename TPatch::ConstIterator  PatchKernelIteratorType;
    typedef typename TPatch::PixelType      KernelPixelType;
    SVMRegressionGrading(const Self&); //purposely not implemented
    
    svm_problem   svm_prob;
    svm_node    * svm_x_space;
    double      * svm_prob_estimates;
    int         * svm_labels;
    svm_node    * svm_x;
    struct svm_model   * svm_model;
    int           svm_j,svm_i;
    
    SVMRegressionGrading() 
    {
      svm_prob.y  = NULL;
      svm_prob.x  = NULL;
      svm_x_space = NULL;
      svm_prob_estimates=NULL;
      svm_labels = NULL;
      svm_x      = NULL;
      svm_model  = NULL;
    
      svm_param.svm_type = C_SVC;
      svm_param.kernel_type = LINEAR;

      svm_param.degree     = 3;
      svm_param.gamma      = 1.0;
      svm_param.coef0      = 0;
      svm_param.nu         = 0.5;
      svm_param.cache_size = 1000;
      svm_param.C          = 1;
      svm_param.eps        = 1e-3;
      svm_param.p          = 0.1;
      svm_param.shrinking  = 1;
      svm_param.probability = 0; 
      svm_param.nr_weight  = 0;
      
      svm_param.weight_label = NULL;
      svm_param.weight = NULL;
      
      norm_min=0.0;
      norm_max=1.0;
    }
    
    ~SVMRegressionGrading()
    {
      cleanup();
    }
    
    void cleanup(void)
    {
      if(svm_prob.y)
        delete [] svm_prob.y;
      
      if(svm_prob.x)
        delete [] svm_prob.x;
      
      if(svm_x_space)
        delete [] svm_x_space;
      
      if(svm_labels)
        delete [] svm_labels;
      
      if(svm_param.weight)
          svm_destroy_param(&svm_param);
      
      if(svm_model)
        svm_free_and_destroy_model(&svm_model);
      
      svm_prob.y=NULL;
      svm_prob.x=NULL;
      svm_x_space=NULL;
      svm_prob_estimates=NULL;
      svm_labels=NULL;
      svm_x=NULL;
      svm_model=NULL;
      
      svm_param.weight_label = NULL;
      svm_param.weight = NULL;
      
    }

  public:
    Pointer clone(void) const
    {
      Pointer _clone=Self::New();
      //TODO: finish 
      //_clone->...=...
      
      _clone->svm_param=svm_param;
      _clone->norm_min=norm_min;
      _clone->norm_max=norm_max;
      
      return _clone;
    }
    
    void allocate_memory(size_t sample_size, size_t feature_size,size_t label_cnt)
    {
      cleanup();
      svm_param.gamma = 1.0/(double)feature_size;
      
      svm_prob.y  =     new double   [ sample_size ];
      svm_prob.x  =     new svm_node*[ sample_size ];
      svm_x_space =     new svm_node [ sample_size * (feature_size + 1) ];
      svm_x       =     new svm_node [ feature_size + 1 ];

      svm_prob.l      = sample_size;

      // pre-allocate features
      for(int i=0,j=0; i<sample_size; i++)
      {
        svm_prob.x[i] = &svm_x_space[j];

        for(int k=0; k<feature_size; k++)
        {
          svm_x_space[j++].index=k+1;
        }
        svm_x_space[j++].index=-1;
      }
      
      for(int k=0; k<feature_size; k++)
      {
        svm_x[k].index = k+1;
        svm_x[k].value = 0.0;
      }
      
      svm_x[feature_size].index=-1;
      
      svm_prob_estimates = new double[label_cnt];
      svm_labels         = new    int[label_cnt];
    }
    
    
    void  start_new_regression(void)
    {
      //TODO: reset internal buffers?
      if(svm_param.weight)
          svm_destroy_param(&svm_param);
      
      if(svm_model)
        svm_free_and_destroy_model(&svm_model);
      
      svm_model=NULL;
      
      svm_param.weight_label = NULL;
      svm_param.weight = NULL;
      
      svm_j=0;
      svm_i=0;
    }
    
    
    void set_input(const NeighborhoodIteratorType &patchIt1,
                   const PatchKernelIteratorType patchKernelBegin,
                   const PatchKernelIteratorType patchKernelEnd
                  )
    {
      
      double norm_k=2.0/(norm_max-norm_min);
      ::size_t i,j=0;
      PatchKernelIteratorType kernel_it;
      //::size_t center = (::size_t) (patchIt1.Size() / 2); // get offset of center pixel
      
      for( i=0, kernel_it=patchKernelBegin; kernel_it<patchKernelEnd; ++kernel_it, ++i )
      {
        // if structuring element is positive, use the pixel under that element
        // in the image
        if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero  )
        {
          svm_x[j++].value=(patchIt1.GetPixel(i)-norm_min)*norm_k-1.0;
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
      double norm_k=2.0/(norm_max-norm_min);
      ::size_t i;
      PatchKernelIteratorType kernel_it;
      svm_prob.y[i]=label;
      
      for( i=0, kernel_it=patchKernelBegin; kernel_it<patchKernelEnd; ++kernel_it, ++i )
      {
        // if structuring element is positive, use the pixel under that element
        // in the image
        
        if( *kernel_it > itk::NumericTraits<KernelPixelType>::Zero  )
        {
          svm_x_space[svm_j].value=(patchIt2.GetPixel(i)-norm_min)*norm_k-1.0;
          svm_j++;
        }
      }
      svm_i++;
    }

    
    bool regress(int &label,double &confidence,std::vector<double> &prob,double &grading)
    {
      svm_model=svm_train(&svm_prob,&svm_param);
      
      if(  svm_param.probability )
      {
        int nr_class=svm_get_nr_class(svm_model);
        svm_get_labels(svm_model,svm_labels);
        
        label = (int)svm_predict_probability(svm_model,svm_x,svm_prob_estimates);
        
        for(int k=0;k<nr_class;k++)
        {
          prob[svm_labels[k]]=svm_prob_estimates[k];
          
          if(svm_labels[k]==label )
            confidence=svm_prob_estimates[k];
        }
      } else {
        label = (int)svm_predict(svm_model,svm_x);
        prob[label] = 1.0;
        confidence = 1.0;
      }
      grading = svm_get_nr_class(svm_model);
      
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

#endif //_ITK_SVM_REGRESSION_GRADING_h_