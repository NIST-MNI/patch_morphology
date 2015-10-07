#ifndef __mincMinimalDistanceNonLocalFilter_h
#define __mincMinimalDistanceNonLocalFilter_h

#include "itkNonLocalPatchesFilter.h"
#include "itkPatchDistance.h"
#include "itkPatchCostFunction.h"

namespace itk
{


  template<class TInputImage, class TOutputImage, class TSearch, class TPatch,class TDistance=L2PatchDistance<TInputImage,TPatch>,class TPreselectionFilter=NOOPPreselection<3> >
  class  MinimalDistanceNonLocalFilter : 
      public NonLocalPatchesFilter<TInputImage, TOutputImage, TSearch, TPatch, TPreselectionFilter>
  {
  public:
    /** Standard class typedefs. */
    typedef MinimalDistanceNonLocalFilter Self;
    typedef NonLocalPatchesFilter<TInputImage, TOutputImage, TSearch, TPatch, TPreselectionFilter>
                                      Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;
    
    /** Standard New method. */
    itkNewMacro(Self);  

    /** Runtime information support. */
    itkTypeMacro(MinimalDistanceNonLocalFilter, 
                NonLocalPatchesFilter);
    
    /** Declaration of pixel type. */
    typedef typename Superclass::PixelType PixelType;

    /** Kernel (structuring element) iterator. */
    typedef typename Superclass::SearchKernelIteratorType SearchKernelIteratorType;
    
    /** Kernel (structuring element) iterator. */
    typedef typename Superclass::PatchKernelIteratorType  PatchKernelIteratorType;

    /** Neighborhood iterator type. */
    typedef typename Superclass::NeighborhoodIteratorType NeighborhoodIteratorType;

    /** Kernel typedef. */
    typedef typename Superclass::SearchKernelType SearchKernelType;
    typedef typename Superclass::PatchKernelType  PatchKernelType;
    
    typedef typename Superclass::PreselectionFilterPointerType PreselectionFilterPointerType;
    
    typedef typename TDistance::Pointer DistancePointerType;

    /** Default boundary condition type */
    typedef typename Superclass::DefaultBoundaryConditionType DefaultBoundaryConditionType;

    /** ImageDimension constants */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TOutputImage::ImageDimension);
    itkStaticConstMacro(KernelDimension, unsigned int,
                        TSearch::NeighborhoodDimension);

    /** Type of the pixels in the Kernel. */
    typedef typename TSearch::PixelType            KernelPixelType;


    /** Set distance functor. */
    itkSetMacro(Distance, DistancePointerType);
    /** Get distance functor. */
    itkGetMacro(Distance, DistancePointerType);


  // #ifdef ITK_USE_CONCEPT_CHECKING
  //   /** Begin concept checking */
  //   itkConceptMacro(InputConvertibleToOutputCheck,
  //     (Concept::Convertible<PixelType, typename TOutputImage::PixelType>));
  //   itkConceptMacro(SameDimensionCheck1,
  //      (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
  //   itkConceptMacro(SameDimensionCheck2,
  //     (Concept::SameDimension<InputImageDimension, KernelDimension>));
  //   itkConceptMacro(InputGreaterThanComparableCheck,
  //     (Concept::GreaterThanComparable<PixelType>));
  //   itkConceptMacro(KernelGreaterThanComparableCheck,
  //     (Concept::GreaterThanComparable<KernelPixelType>));
  //   /** End concept checking */
  // #endif

  protected:
    MinimalDistanceNonLocalFilter();
    ~MinimalDistanceNonLocalFilter() {};

    /** Evaluate image neighborhood with kernel to find the new value 
    * for the center pixel value
    *
    */

    PixelType Evaluate(const NeighborhoodIteratorType &searchIt,
                      const  NeighborhoodIteratorType &patchIt1,
                      NeighborhoodIteratorType &patchIt2,
                      const SearchKernelIteratorType searchKernelBegin,
                      const SearchKernelIteratorType searchKernelEnd,
                      const PatchKernelIteratorType patchKernelBegin,
                      const PatchKernelIteratorType patchKernelEnd,
                      PreselectionFilterPointerType flt);


    void PrintSelf(std::ostream& os, itk::Indent indent) const;

  private:
    MinimalDistanceNonLocalFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
    
    DistancePointerType m_Distance;
  }; // end of class
  
} //minc


#include "itkMinimalDistanceNonLocalFilter.txx"

#endif //__mincMinimalDistanceNonLocalFilter_h

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;show-tabs on;tab-width 2;hl c++
