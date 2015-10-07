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
#ifndef __mincVariableNoiseNonLocalFilter_h
#define __mincVariableNoiseNonLocalFilter_h

#include "itkNonLocalPatchesFilter.h"
#include "itkPatchDistance.h"
#include "itkPatchCostFunction.h"

namespace itk
{
  template<class TInputImage, 
          class TOutputImage, 
          class TSearch, 
          class TPatch,
          class TDistance=L2PatchDistance<TInputImage,TPatch>,
          class TWeight=InvExpWeight<double>,
          class TPreselectionFilter=NOOPPreselection<3> >
  class  VariableNoiseNonLocalFilter : 
      public itk::ImageToImageFilter<TInputImage, TOutputImage>
  {
  public:
    /** Standard class typedefs. */
    typedef VariableNoiseNonLocalFilter Self;
    typedef itk::ImageToImageFilter<TInputImage, TOutputImage>  Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;
    
    /** Standard New method. */
    itkNewMacro(Self);  

    /** Runtime information support. */
    itkTypeMacro(VariableNoiseNonLocalFilter,itk::ImageToImageFilter);
    
    /** Image related typedefs. */
    typedef TInputImage                                   InputImageType;
    typedef TOutputImage                                  OutputImageType;
    typedef typename TInputImage::RegionType              RegionType;
    typedef typename TInputImage::SizeType                SizeType;
    typedef typename TInputImage::IndexType               IndexType;
    typedef typename TInputImage::PixelType               PixelType;
    typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;
    typedef typename TOutputImage::Pointer            		OutputImagePointer;
    
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
    typedef itk::ConstNeighborhoodIterator<TInputImage>  	NeighborhoodIteratorType;

    /** Kernel typedef. */
    typedef TSearch SearchKernelType;
    typedef TPatch  PatchKernelType;

    typedef typename TPreselectionFilter::Pointer PreselectionFilterPointerType;
    typedef typename TDistance::Pointer           DistancePointerType;

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
    
    /** VariableNoiseNonLocalFilter need to make sure they request enough of an
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
    itkSetMacro(Distance, DistancePointerType);
    /** Get distance functor. */
    itkGetMacro(Distance, DistancePointerType);

    /** Set Input Image */
    void SetInput1( const TInputImage * image1);

    /** Set Noise Image */
    void SetInput2( const TInputImage * image2);

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
    
    /** assign ROI Image */
    itkSetMacro(RoiImage,RoiImagePointer);
    /** get ROI Image */
    itkGetMacro(RoiImage,RoiImagePointer);
    
    /** assign preselection filter */
    itkSetMacro(PreselectionFilter,PreselectionFilterPointerType);
    /** get preselection filter */
    itkGetMacro(PreselectionFilter,PreselectionFilterPointerType);

    /** set output weight flag */
    itkSetMacro(OutputMeanWeight,bool);
    /** get output weight flag */
    itkGetMacro(OutputMeanWeight,bool);
    
    /** set beta parameter (abjust smoothing) */
    itkSetMacro(Beta,double);
    /** get beta parameter (abjust smoothing) */
    itkGetMacro(Beta,double);
    
  protected:
    VariableNoiseNonLocalFilter();
    ~VariableNoiseNonLocalFilter() {};

    /** prepare before starting the thread */
    void BeforeThreadedGenerateData(void);

    /** Multi-thread version GenerateData. */
  #if (ITK_VERSION_MAJOR==3)
      void  ThreadedGenerateData (const OutputImageRegionType& outputRegionForThread,
                                  int threadId);
  #else
      void  ThreadedGenerateData (const OutputImageRegionType& outputRegionForThread,
                                  itk::ThreadIdType threadId);
  #endif    
    
    
    /** Evaluate image neighborhood with kernel to find the new value 
    * for the center pixel value
    *
    */
    
    PixelType Evaluate(const NeighborhoodIteratorType &searchIt,
                      double sigma,
                      const  NeighborhoodIteratorType &patchIt1,
                      NeighborhoodIteratorType &patchIt2,
                      const SearchKernelIteratorType searchKernelBegin,
                      const SearchKernelIteratorType searchKernelEnd,
                      const PatchKernelIteratorType patchKernelBegin,
                      const PatchKernelIteratorType patchKernelEnd,
                      PreselectionFilterPointerType flt
                      );

    
    void PrintSelf(std::ostream& os, itk::Indent indent) const;

  private:
    VariableNoiseNonLocalFilter(const Self&); //purposely not implemented
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
    
    DistancePointerType m_Distance;
    TWeight   m_Weight;
    bool      m_OutputMeanWeight;
    RoiImagePointer    m_RoiImage;
    PreselectionFilterPointerType m_PreselectionFilter;
    double    m_Beta;

  }; // end of class

} //minc

#include "itkVariableNoiseNonLocalFilter.txx"

#endif //__mincVariableNoiseNonLocalFilter_h

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl C++
