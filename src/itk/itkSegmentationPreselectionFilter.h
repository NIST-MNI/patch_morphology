#ifndef __mincSegmentationPreselectionFilter_h__
#define __mincSegmentationPreselectionFilter_h__

#include <vector>

#include <itkLightObject.h>

#include "itkLocalMeanFilter.h"
#include "itkLocalSDFilter.h"
#include <itkProgressAccumulator.h>
#include <itkImageDuplicator.h>

#include <itkLightObject.h>

namespace itk
{

  //! preselect all filter
  template <int VImageDimension=3> class NOOPSegmentationPreselection: public itk::LightObject
  {
    typedef typename itk::ImageBase< VImageDimension > ImageBase;
    typedef typename ImageBase::IndexType IndexType;
  public:
    typedef NOOPSegmentationPreselection             Self;
    itkTypeMacro(NOOPSegmentationPreselection, itk::LightObject);    

    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    itkNewMacro(Self);
    
    bool select(const IndexType& src,size_t trg_idx,const IndexType& trg)
    {
      return true;
    }
    
    bool operator!=(const NOOPSegmentationPreselection &a)
    {
      return false;
    }
    
    NOOPSegmentationPreselection& operator=(const NOOPSegmentationPreselection& another)
    {
      return *this;
    }
    
    void SetThreshold(double b)
    {
      //NOOP
    }
    
    double GetThreshold(void) 
    {
      return 0.0;
    }
    
  };

  template<int dim>
      std::ostream & operator<<(std::ostream &os, const NOOPSegmentationPreselection<dim> &w)
  {
    os << "[ NOOPSegmentationPreselection<dim> ]";
    return os;
  }

  template <class TImage,class TPatch> class SegmentationMeanAndSdPreselectionFilter: public itk::LightObject
  {
  public:
    typedef SegmentationMeanAndSdPreselectionFilter Self;
    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    /** Standard New method. */
    itkNewMacro(Self);
    /** Runtime information support. */
    itkTypeMacro(SegmentationMeanAndSdPreselectionFilter,itk::LightObject);
    
    typedef TImage               ImageType;
    typedef typename ImageType::IndexType IndexType;
    typedef typename ImageType::Pointer   ImagePointer;
    typedef typename ImageType::ConstPointer   ConstImagePointer;
    typedef typename ImageType::PixelType  PixelType;
    typedef std::vector<ImagePointer> TargetImageVector;
    
    typedef LocalSDFilter<TImage,TImage,TPatch>   SDFilterType;
    typedef LocalMeanFilter<TImage,TImage,TPatch> MeanFilterType;	
    
    typedef typename SDFilterType::Pointer 		    SDFilterPointer;
    typedef typename MeanFilterType::Pointer 		  MeanFilterPointer;
    
    typedef itk::ImageDuplicator<TImage>          DuplicatorType;
    typedef typename DuplicatorType::Pointer      DuplicatorPointer;
    
  public:
    
    //! evaluate if given indexes are similar enough, according to Pierrick's code
    bool select(const IndexType& src,size_t trg_idx,const IndexType& trg)
    {
      const double epsi=1e-4;
      double Mean=m_src_mean->GetPixel(src);
      double Var=m_src_sd->GetPixel(src);
      double TMean=m_target_mean[trg_idx]->GetPixel(trg);
      double TVar=m_target_sd[trg_idx]->GetPixel(trg);
      
      /*Similar Luminance and contrast -> Cf Wang TIP 2004*/
      double th = ((2 * Mean * TMean + epsi) / ( Mean*Mean + TMean*TMean + epsi))  * ((2 * Var*TVar + epsi) / (Var*Var + TVar*TVar + epsi));
      
      return th>m_Threshold;
    }
    
    bool operator!=(const SegmentationMeanAndSdPreselectionFilter &a)
    {
      return true;
    }
    
    SegmentationMeanAndSdPreselectionFilter& operator=(const SegmentationMeanAndSdPreselectionFilter& another)
    {
      m_Patch=another.m_Patch;
      
      m_src_mean=another.m_src_mean;
      m_src_sd=another.m_src_sd;
      
      m_target_mean=another.m_target_mean;
      m_target_sd=another.m_target_sd;
      
      m_Threshold=another.m_Threshold;
      
      return *this;
    }
    
    SegmentationMeanAndSdPreselectionFilter():
      m_Threshold(0.97)
    {
    }
    
    SegmentationMeanAndSdPreselectionFilter(const SegmentationMeanAndSdPreselectionFilter& another):
      m_Patch(another.m_Patch),
      m_src_mean(another.m_src_mean),
      m_src_sd(another.m_src_sd),
      m_target_mean(another.m_target_mean),
      m_target_sd(another.m_target_sd),
      m_Threshold(another.m_Threshold)
    {
    }
    
    void SetPatchKernel(const TPatch& p)
    {
      m_Patch=p;
    }
    
    TPatch& GetPatchKernel(void)
    {
      return m_Patch;
    }
    
    void SetThreshold(PixelType b)
    {
      m_Threshold=b;
    }
    
    PixelType GetThreshold(void) 
    {
      return m_Threshold;
    }
    
    
    void SetImages(ImagePointer src,TargetImageVector& trg)
    {
      
      //itk::ProgressAccumulator::Pointer progress = itk::ProgressAccumulator::New();
      //progress->SetMiniPipelineFilter(this);	
      
      SDFilterPointer   sd_filter=SDFilterType::New();
      MeanFilterPointer mean_filter=MeanFilterType::New();
      DuplicatorPointer dup=DuplicatorType::New();
      
      sd_filter->SetKernel( m_Patch );
      mean_filter->SetKernel( m_Patch );
      
      m_target_mean.clear();
      m_target_sd.clear();
      
      //progress->RegisterInternalFilter(m_target_mean,0.5/(trg.size()+1));
      //progress->RegisterInternalFilter(m_target_sd,0.5/(trg.size()+1));
      
      for(size_t i=0;i<trg.size();i++)
      {

        mean_filter->SetInput(trg[i]);
        sd_filter->SetInput(trg[i]);
        
        mean_filter->Update();
        sd_filter->Update();
        dup->SetInputImage(mean_filter->GetOutput());dup->Update();
        m_target_mean.push_back(dup->GetOutput());
        
        dup->SetInputImage(sd_filter->GetOutput());dup->Update();
        m_target_sd.push_back(dup->GetOutput());
      }
      
      mean_filter->SetInput(src);
      sd_filter->SetInput(src);
      
      mean_filter->Update();
      sd_filter->Update();	
      
      dup->SetInputImage(mean_filter->GetOutput());dup->Update();
      m_src_mean=dup->GetOutput();
      
      dup->SetInputImage(sd_filter->GetOutput());dup->Update();
      m_src_sd=dup->GetOutput();
    }
    
  protected:
    TPatch    m_Patch;
    PixelType m_Threshold;
    
    TargetImageVector m_target_mean;
    TargetImageVector m_target_sd;
    
    ImagePointer      m_src_mean;
    ImagePointer      m_src_sd;
  };

  template<class TImage,class TPatch>
      std::ostream & operator<<(std::ostream &os, const SegmentationMeanAndSdPreselectionFilter<TImage,TPatch> &w)
  {
    os << "[ SegmentationMeanAndSdPreselectionFilter<TImage,TPatch> ]";
    return os;
  }
}; //minc

#endif //__mincSegmentationPreselectionFilter_h__
