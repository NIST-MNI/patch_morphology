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
#ifndef __itkPreselectionFilter_h__

#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkLightObject.h>

#include "itkLocalMeanFilter.h"
#include "itkLocalSDFilter.h"


namespace itk
{
  //! preselect all filter
  template <int VImageDimension=3> class NOOPPreselection: public itk::LightObject
  {
    typedef typename itk::ImageBase< VImageDimension > ImageBase;
    typedef typename ImageBase::IndexType IndexType;
  public:
    typedef NOOPPreselection             Self;
    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    itkNewMacro(Self);
    itkTypeMacro(NOOPPreselection, itk::LightObject);    

    bool select(const IndexType& src,const IndexType& trg)  const
    {
      return true;
    }
    
    bool operator!=(const NOOPPreselection &a)  const
    {
      return false;
    }
    
    NOOPPreselection& operator=(const NOOPPreselection& another)
    {
      return *this;
    }

    void SetImage(DataObject::Pointer src)
    {
      //NOOP
    }
    
     void SetThresholds(float th_m,float th_v)
     {
      //NOOP
     }
  };
  
  template<int dim>
      std::ostream & operator<<(std::ostream &os, const NOOPPreselection<dim> &w)
  {
    os << "[ NOOPPreselection<dim> ]";
    return os;
  }
  
  //! random preselection filter
  template <int VImageDimension=3> class RandomPreselection: public itk::LightObject
  {
    typedef typename itk::ImageBase< VImageDimension > ImageBase;
    typedef typename ImageBase::IndexType IndexType;
  protected:
    unsigned long _probability;
    itk::Statistics::MersenneTwisterRandomVariateGenerator _rng;
  public:

    typedef RandomPreselection Self;
    
    itkTypeMacro(RandomPreselection, itk::LightObject);    
    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    itkNewMacro(Self);
    
    
    RandomPreselection(double probability=0.5)
    {
      _probability=probability*8589934591;//for integer comparision
    }
    
    bool select(const IndexType& src,const IndexType& trg)
    {
      return _rng.GetIntegerVariate()<_probability;
    }
    
    bool operator!=(const RandomPreselection &a)  const
    {
      return false;
    }
    
    RandomPreselection& operator=(const RandomPreselection& another)
    {
      _probability=another._probability;
      return *this;
    }
    
    void SetImage(DataObject::Pointer src)
    {
      //NOOP
    }
    
    void SetThresholds(float th_m,float th_v)
    {
    //NOOP
    }
    
  };
  
  
  template <class TImage,class TPatch> class MeanAndSdPreselectionFilter: public itk::LightObject
  {
  public:
    typedef MeanAndSdPreselectionFilter Self;
    
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    /** Standard New method. */
    itkNewMacro(Self);
    /** Runtime information support. */
    itkTypeMacro(MeanAndSdPreselectionFilter,itk::LightObject);
    
    typedef TImage               ImageType;
    typedef typename ImageType::IndexType IndexType;
    typedef typename ImageType::Pointer   ImagePointer;
    typedef typename ImageType::ConstPointer   ConstImagePointer;
    typedef typename ImageType::PixelType  PixelType;
    typedef std::vector<ImagePointer> TargetImageVector;
    
    typedef LocalSDFilter<TImage,TImage,TPatch>   SDFilterType;
    typedef LocalMeanFilter<TImage,TImage,TPatch> MeanFilterType; 
    
    typedef typename SDFilterType::Pointer        SDFilterPointer;
    typedef typename MeanFilterType::Pointer      MeanFilterPointer;
    
//     typedef itk::ImageDuplicator<TImage>          DuplicatorType;
//     typedef typename DuplicatorType::Pointer      DuplicatorPointer;
    
  public:
    
    //! evaluate if given indexes are similar enough, according to Pierrick's code
    bool select(const IndexType& src,const IndexType& trg)  const
    {
      const double epsi=1e-4;
      
      double SMean= this->m_src_mean->GetPixel(src);
      double SVar=  this->m_src_sd->GetPixel(src);
      
      double TMean= this->m_src_mean->GetPixel(trg);
      double TVar=  this->m_src_sd->GetPixel(trg);
      
      /*Similar Luminance and contrast -> Cf Wang TIP 2004*/
//       double th = ((2 * Mean * TMean + epsi) / ( Mean*Mean + TMean*TMean + epsi))  * ((2 * Var*TVar + epsi) / (Var*Var + TVar*TVar + epsi));
      
//       return th>m_Threshold;

      if(TMean < epsi || TVar<epsi) return true;

      double ratio  = SMean/TMean;
      double ratio2 = SVar/TVar;
      
      return  (m_th_m <= ratio)  &&  (ratio <= m_th_m_inv) && 
              (m_th_v <= ratio2) && (ratio2 <= m_th_v_inv) ;
    }
    
    bool operator!=(const MeanAndSdPreselectionFilter &a)  const
    {
      return true;
    }
    
    MeanAndSdPreselectionFilter& operator=(const MeanAndSdPreselectionFilter& another)
    {
      this->m_Patch=another.m_Patch;
      
      this->m_src_mean=another.m_src_mean;
      this->m_src_sd=another.m_src_sd;
      
      this->SetThresholds(another.m_th_m,another.m_th_v);
      return *this;
    }
    
    MeanAndSdPreselectionFilter()
    {
      this->SetThresholds(0.95,0.5);
    }
    
    MeanAndSdPreselectionFilter(const MeanAndSdPreselectionFilter& another):
      m_Patch(another.m_Patch),
      m_src_mean(another.m_src_mean),
      m_src_sd(another.m_src_sd)
    {
      this->SetThresholds(another.m_th_m,another.m_th_v);
    }
    
    void SetPatchKernel(const TPatch& p)
    {
      m_Patch=p;
    }
    
    TPatch& GetPatchKernel(void)
    {
      return m_Patch;
    }
    
    const TPatch& GetPatchKernel(void) const
    {
      return m_Patch;
    }
    
    void SetThresholds(PixelType th_m,PixelType th_v)
    {
      m_th_m=th_m;
      m_th_v=th_v;
      
      m_th_m_inv=1.0/m_th_m;
      m_th_v_inv=1.0/m_th_v;
    }
    
    void SetImage(ImagePointer src)
    {
      SDFilterPointer   sd_filter   = SDFilterType::New();
      MeanFilterPointer mean_filter = MeanFilterType::New();
      
      sd_filter->SetKernel( m_Patch );
      mean_filter->SetKernel( m_Patch );
      
      mean_filter->SetInput(src);
      sd_filter->SetInput(src);
      
      mean_filter->Update();
      sd_filter->Update();  
      
      m_src_mean=mean_filter->GetOutput();
      m_src_sd=sd_filter->GetOutput();
      
      m_src_mean->DisconnectPipeline  (   ) ;
      m_src_sd->DisconnectPipeline  (   ) ;
    }
    
    PixelType Get_th_m(void) const
    {
      return m_th_m;
    }
    
    PixelType Get_th_v(void) const
    {
      return m_th_v;
    }
    
  protected:
    TPatch            m_Patch;
    PixelType         m_th_m,m_th_m_inv;
    PixelType         m_th_v,m_th_v_inv;

    ImagePointer      m_src_mean;
    ImagePointer      m_src_sd;
  };

  template<class TImage,class TPatch>
      std::ostream & operator<<(std::ostream &os, const MeanAndSdPreselectionFilter<TImage,TPatch> &w)
  {
    os << "[ MeanAndSdPreselectionFilter<TImage,TPatch> ";
    os << "m_th_m="<<w.Get_th_m()<<" ";
    os << "m_th_v="<<w.Get_th_v()<<" ";
    os << "Patch="<<w.GetPatchKernel() << " ]"<<std::endl;
    return os;
  }  
};

#endif //__itkPreselectionFilter_h__

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl C++
