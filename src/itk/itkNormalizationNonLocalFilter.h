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
#ifndef __mincNormalizationNonLocalFilter_h
#define __mincNormalizationNonLocalFilter_h

#include <itkLightObject.h>
#include "itkNonLocalPatchesFilter.h"
#include "itkPatchDistance.h"
#include "itkPatchCostFunction.h"
#include "itkNormalizationProcess.h"
#include "itkSegmentationPreselectionFilter.h"

namespace itk
{
  template<class InputPixel,class OutputPixel>
      std::ostream & operator<<(std::ostream &os, const MedianNormalization<InputPixel,OutputPixel> &w)
  {
    os << "[ MedianNormalization<InputPixel,OutputPixel> ]";
    return os;
  }

  template<class TFeatureImage, 
          class TOutputImage, 
          class TSearch, 
          class TPatch,
          class TSampleProcess,
          class TDistance=L2PatchDistance<TFeatureImage,TPatch>,
          class TPreselectionFilter=NOOPSegmentationPreselection<3>,
          class TRealType=float>
  class  NormalizationNonLocalFilter : 
      public itk::ImageToImageFilter<TFeatureImage, TOutputImage>
  {
    public:
    
    /** Standard class typedefs. */
    typedef NormalizationNonLocalFilter Self;
    typedef itk::ImageToImageFilter<TFeatureImage, TOutputImage>  Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;

    /** Standard New method. */
    itkNewMacro(Self);  

    /** Runtime information support. */
    itkTypeMacro(NormalizationNonLocalFilter,itk::ImageToImageFilter);

    /** Image related typedefs. */
    typedef TFeatureImage                                   FeatureImageType;
    typedef TOutputImage                                    OutputImageType;

    typedef typename TFeatureImage::RegionType              RegionType;
    typedef typename TFeatureImage::SizeType                SizeType;
    typedef typename TFeatureImage::IndexType               IndexType;
    typedef typename TFeatureImage::PixelType               PixelType;
    typedef typename TOutputImage::PixelType                OutputPixelType;
    
    typedef typename Superclass::OutputImageRegionType      OutputImageRegionType;
    typedef typename TOutputImage::Pointer                  OutputImagePointerType;
    
    typedef typename std::vector<typename TFeatureImage::Pointer >            NormalizationLibraryType;
    
    typedef typename itk::Image<unsigned char, TFeatureImage::ImageDimension> RoiImageType;
    
    typedef typename RoiImageType::Pointer                RoiImagePointer;

    
    /** Image related typedefs. */
    itkStaticConstMacro(ImageDimension, unsigned int, TFeatureImage::ImageDimension);
    typedef typename std::vector<NormSampleDistance<PixelType> > NormalizationSampleDistance;

    /** Parameters typedefs */
    typedef typename TSampleProcess::Pointer SampleProcessPointer;
    typedef typename TDistance::Pointer      DistancePointer;
    typedef typename TPreselectionFilter::Pointer PreselectionFilterPointerType;
    
    /** Typedef for boundary conditions. */
    typedef itk::ImageBoundaryCondition<FeatureImageType> *       ImageBoundaryConditionPointerType;
    typedef itk::ImageBoundaryCondition<FeatureImageType> const * ImageBoundaryConditionConstPointerType;
    typedef itk::ConstantBoundaryCondition<FeatureImageType>      DefaultBoundaryConditionType;
    

  /** Neighborhood iterator type. */
    typedef itk::ConstNeighborhoodIterator<TFeatureImage>   NeighborhoodIteratorType;
    
    /** Kernel typedef. */
    typedef TSearch SearchKernelType;
    typedef TPatch  PatchKernelType;

    /** Internal float image type. */
    typedef itk::Image<TRealType,itkGetStaticConstMacro(ImageDimension)> FloatImageType;
    typedef typename FloatImageType::Pointer         FloatImagePointerType;
    
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
                        TFeatureImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TOutputImage::ImageDimension);
    itkStaticConstMacro(KernelDimension, unsigned int,
                        TSearch::NeighborhoodDimension);

    /** Type of the pixels in the Kernel. */
    typedef typename TSearch::PixelType            KernelPixelType;
    
    /** NormalizationNonLocalFilter need to make sure they request enough of an
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
    
    void SetNormalizationLibrary(NormalizationLibraryType &s)
    {
      m_NormalizationLibrary=s;
    }
    
    NormalizationLibraryType & GetNormalizationLibrary(void)
    {
      return m_NormalizationLibrary;
    }

    /** Set distance functor. */
    itkSetMacro(Distance, DistancePointer);
    /** Get distance functor. */
    itkGetMacro(Distance, DistancePointer);
    
    /** assign ROI Image */
    itkSetMacro(RoiImage,RoiImagePointer);
    /** get ROI Image */
    itkGetMacro(RoiImage,RoiImagePointer);
    
    /** assign ROI Image */
    itkSetMacro(Process,SampleProcessPointer);
    /** get ROI Image */
    itkGetMacro(Process,SampleProcessPointer);
    
    /** assign preselection filter */
    itkSetMacro(PreselectionFilter,PreselectionFilterPointerType);
    /** get preselection filter */
    itkGetMacro(PreselectionFilter,PreselectionFilterPointerType);
    
    /** assign confidence threshold */
    itkSetMacro(ConfidenceThreshold,TRealType);
    /** get confidence threshold */
    itkGetMacro(ConfidenceThreshold,TRealType);
    
    /** get number of voxels below confidence in last iteration */
    itkGetMacro(NonConfidentVoxels,TRealType);

    /** assign preload results */
    itkSetMacro(PreloadedResults,OutputImagePointerType);
    
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
    NormalizationNonLocalFilter();
    ~NormalizationNonLocalFilter() {};

    /** prepare before starting the thread */
    void BeforeThreadedGenerateData(void);

    /** Multi-thread version GenerateData. */
    void  ThreadedGenerateData (const OutputImageRegionType& outputRegionForThread,
                                itk::ThreadIdType threadId);
    /** finish calculations */
    void AfterThreadedGenerateData(void);
    
    

    void PrintSelf(std::ostream& os, itk::Indent indent) const;

  private:
    NormalizationNonLocalFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
    
    /** kernel or structuring element to use. */
    SearchKernelType        m_SearchKernel;
    /** kernel or structuring element to use. */
    PatchKernelType         m_PatchKernel;
    /** Pointer to a persistent boundary condition object used
    * for the image iterator. */
    ImageBoundaryConditionPointerType 
                            m_BoundaryCondition;
    /** Default boundary condition */
    DefaultBoundaryConditionType 
                            m_DefaultBoundaryCondition;
    /** Patch distance functor */
    DistancePointer         m_Distance;
    /** sample list processing functor */
    SampleProcessPointer    m_Process;
    /** Normalization library */
    NormalizationLibraryType m_NormalizationLibrary;
    /** ROI mask */
    RoiImagePointer         m_RoiImage;
    /** Preselection filter */
    PreselectionFilterPointerType     m_PreselectionFilter;
    
    TRealType               m_ConfidenceThreshold;
    /** Number of voxels which had confidence below threshold in last iteration */
    TRealType               m_NonConfidentVoxels;
    /** Number of voxels which had confidence below threshold in last iteration (in one thread)*/
    std::vector<TRealType>  m_ThreadNonConfidentVoxels;
    
    /** preloaded results, to initialize segmentation */
    OutputImagePointerType  m_PreloadedResults;
  }; // end of class

} //minc

#include "itkNormalizationNonLocalFilter.txx"

#endif //__mincNormalizationNonLocalFilter_h