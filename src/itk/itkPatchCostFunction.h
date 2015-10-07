#ifndef _mincPatchCostFunctions_h_
#define _mincPatchCostFunctions_h_

#include <iostream>
#include <itkLightObject.h>

namespace itk
{
  template<class T=double> class InvExpWeight
  {
  public:
    double m_Beta;
    
    InvExpWeight():
      m_Beta(1.0)
    {
      
    }
    
    double operator()(const T& distance,const T& sigma2) const
    {
      return ::exp(-distance/(2*sigma2*m_Beta));
    }
    
    bool operator!=(const InvExpWeight &a)
    {
      return false;
    }
    
    void SetBeta(double beta)
    {
      m_Beta=beta;
    }
    
    double GetBeta(void) const
    {
      return m_Beta;
    }
  };

  template<class T>
      std::ostream & operator<<(std::ostream &os, const InvExpWeight<T> &w)
  {
    os << "[ InvExpWeight<T> beta=" << w.GetBeta() <<" ]";
    return os;
  }

  template <class TImage,class TPatch> 
        std::ostream & operator<<(std::ostream &os, const L2PatchDistance<TImage,TPatch> &d)
  {
    os << "[ L2PatchDistance<TImage,TPatch> ]";
    return os;
  }

}//minc


#endif //_mincPatchCostFunctions_h_

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl C++
