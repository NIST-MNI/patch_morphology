#ifndef __mincPreselectionFilter_h__

#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkLightObject.h>


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

    bool select(const IndexType& src,const IndexType& trg) 
    {
      return true;
    }
    
    bool operator!=(const NOOPPreselection &a)
    {
      return false;
    }
    
    NOOPPreselection& operator=(const NOOPPreselection& another)
    {
      return *this;
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
    
    bool operator!=(const RandomPreselection &a)
    {
      return false;
    }
    
    RandomPreselection& operator=(const RandomPreselection& another)
    {
      _probability=another._probability;
      return *this;
    }
  };
};

#endif //__mincPreselectionFilter_h__

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl C++
