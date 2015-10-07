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
#ifndef __mincNoiseMedianDistanceNonLocalFilter_h
#define __mincNoiseMedianDistanceNonLocalFilter_h

#include "itkNonLocalPatchesFilter.h"
#include "itkPatchDistance.h"
#include "itkPatchCostFunction.h"
#include "itkLocalMeanFilter.h"
#include "itkMedianDistanceNonLocalFilter.h"
#include <itkSubtractImageFilter.h>
#include <itkDivideImageFilter.h>

namespace itk
{

  template<class TInputImage, class TOutputImage, class TSearch, class TPatch,
          class TDistance=L2PatchDistance<TInputImage,TPatch>,class TWeight=InvExpWeight<double>,
          class TPreselectionFilter=NOOPPreselection<3> >
  class  NoiseMedianDistanceNonLocalFilter : 
      public itk::ImageToImageFilter<TInputImage, TOutputImage>
  {
  public:
    /** Standard class typedefs. */
    typedef NoiseMedianDistanceNonLocalFilter Self;
    typedef itk::ImageToImageFilter<TInputImage, TOutputImage>  Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;
    
    /** Standard New method. */
    itkNewMacro(Self);  

    /** Runtime information support. */
    itkTypeMacro(NoiseMedianDistanceNonLocalFilter,itk::ImageToImageFilter);
    
    /** Image related typedefs. */
    typedef TInputImage                                   InputImageType;
    typedef TOutputImage                                  OutputImageType;
    typedef typename TInputImage::RegionType              RegionType;
    typedef typename TInputImage::SizeType                SizeType;
    typedef typename TInputImage::IndexType               IndexType;
    typedef typename TInputImage::PixelType               PixelType;
    typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;
    typedef typename TOutputImage::Pointer                OutputImagePointer;
    
    /** Image related typedefs. */
    itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

    /** Typedef for boundary conditions. */
    typedef itk::ImageBoundaryCondition<InputImageType> *       ImageBoundaryConditionPointerType;
    typedef itk::ImageBoundaryCondition<InputImageType> const * ImageBoundaryConditionConstPointerType;
    typedef itk::ConstantBoundaryCondition<InputImageType>      DefaultBoundaryConditionType;

  #if (ITK_VERSION_MAJOR==3)  
    typedef typename itk::OrientedImage<unsigned char, TInputImage::ImageDimension> RoiImageType;
  #else
    typedef typename itk::Image<unsigned char, TInputImage::ImageDimension> RoiImageType;
  #endif  
    typedef typename RoiImageType::Pointer                RoiImagePointer;
    

  /** Neighborhood iterator type. */
    typedef itk::ConstNeighborhoodIterator<TInputImage>   NeighborhoodIteratorType;

    /** Kernel typedef. */
    typedef TSearch SearchKernelType;
    typedef TPatch  PatchKernelType;
    
    /** Kernel (structuring element) iterator. */
    typedef typename SearchKernelType::ConstIterator SearchKernelIteratorType;
    
    /** Kernel (structuring element) iterator. */
    typedef typename PatchKernelType::ConstIterator  PatchKernelIteratorType;
    
    /** Set kernel (structuring element). */
    itkSetMacro(SearchKernel, SearchKernelType);
    
    /** Set kernel (structuring element). */
    itkSetMacro(PatchKernel, PatchKernelType);

    /** Get the kernel (structuring element). */
    itkGetConstReferenceMacro(SearchKernel, SearchKernelType);
    
    /** Get the kernel (structuring element). */
    itkGetConstReferenceMacro(PatchKernel, PatchKernelType);

    /** ImageDimension constants */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TOutputImage::ImageDimension);
    itkStaticConstMacro(KernelDimension, unsigned int,
                        TSearch::NeighborhoodDimension);

    /** Type of the pixels in the Kernel. */
    typedef typename TSearch::PixelType            KernelPixelType;
    
    /** NoiseMedianDistanceNonLocalFilter need to make sure they request enough of an
    * input image to account for the structuring element size.  The input
    * requested region is expanded by the radius of the Search+Patch structuring element.
    * If the request extends past the LargestPossibleRegion for the input,
    * the request is cropped by the LargestPossibleRegion. */
    void GenerateInputRequestedRegion();
    
    /** Allows a user to override the internal boundary condition. Care should be
    * be taken to ensure that the overriding boundary condition is a persistent
    * object during the time it is referenced.  The overriding condition
    * can be of a different type than the default type as long as it is
    * a subclass of ImageBoundaryCondition. */
    void OverrideBoundaryCondition(const ImageBoundaryConditionPointerType i)
      { m_BoundaryCondition = i; }

    /** Rest the boundary condition to the default */
    void ResetBoundaryCondition()
      { m_BoundaryCondition = &m_DefaultBoundaryCondition; }
    

    /** Set weight functor. */
    itkSetMacro(Weight, TWeight);
    /** Get weight functor. */
    itkGetMacro(Weight, TWeight);

    /** Set distance functor. */
    itkSetMacro(Distance, TDistance);
    /** Get distance functor. */
    itkGetMacro(Distance, TDistance);

    /** assign ROI Image */
    itkSetMacro(RoiImage,RoiImagePointer);
    /** get ROI Image */
    itkGetMacro(RoiImage,RoiImagePointer);

    
    /** assign preselection filter */
    itkSetMacro(PreselectionFilter,TPreselectionFilter);
    /** get preselection filter */
    itkGetMacro(PreselectionFilter,TPreselectionFilter);
    
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
    NoiseMedianDistanceNonLocalFilter();
    ~NoiseMedianDistanceNonLocalFilter() {};
    
    void GenerateData();
    

    
    void PrintSelf(std::ostream& os, itk::Indent indent) const;

  private:
    NoiseMedianDistanceNonLocalFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
    
    /** kernel or structuring element to use. */
    SearchKernelType m_SearchKernel;
    
    /** kernel or structuring element to use. */
    PatchKernelType m_PatchKernel;

    /** Pointer to a persistent boundary condition object used
    * for the image iterator. */
    ImageBoundaryConditionPointerType m_BoundaryCondition;

    /** Default boundary condition */
    DefaultBoundaryConditionType m_DefaultBoundaryCondition;
    
    TDistance m_Distance;
    TWeight   m_Weight;
    RoiImagePointer m_RoiImage;
    TPreselectionFilter m_PreselectionFilter;


  }; // end of class

} //minc

#include "itkNoiseMedianDistanceNonLocalFilter.txx"

#endif //__mincNoiseMedianDistanceNonLocalFilter_h