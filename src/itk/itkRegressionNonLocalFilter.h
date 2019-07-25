#ifndef __mincRegressionNonLocalFilter_h
#define __mincRegressionNonLocalFilter_h

#include <itkLightObject.h>
#include "itkNonLocalPatchesFilter.h"
#include "itkPatchDistance.h"
#include "itkPatchCostFunction.h"
#include "itkSegmentationProcess.h"
#include "itkSegmentationPreselectionFilter.h"


namespace itk
{

  template<class TFeatureImage, class TLabelImage> struct ImagePair
  {
  public:
    typedef typename TFeatureImage::Pointer FeatureImagePointer;
    typedef typename TLabelImage::Pointer LabelImagePointer;
    
    /** Feature volume */
    FeatureImagePointer feature; 
    /** Label volume*/
    LabelImagePointer   label;
    /** Subject grading value*/
    double              grading;
    /** Subject global similarity*/
    double              similarity;
    /** Subject group */
    int  group;
    
    ImagePair(FeatureImagePointer _f,LabelImagePointer _l,double _g=1.0,double _s=0.0,int _group=-1):
    feature(_f),label(_l),grading(_g),similarity(_s),group(_group)
    {
    }
  };

  template<class TFeatureImage, class TLabelImage> 
  bool operator==(const ImagePair<TFeatureImage,TLabelImage>& a,const ImagePair<TFeatureImage,TLabelImage>& b)
  {
    return false;
  }

  template<class TFeatureImage, class TLabelImage> 
  bool operator!=(const ImagePair<TFeatureImage,TLabelImage>& a,const ImagePair<TFeatureImage,TLabelImage>& b)
  {
    return true;
  }


template<class TFeatureImage, 
        class TLabelImage, 
        class TOutputImage, 
        class TSearch, 
        class TPatch,
        class TRegressor,
        class TPreselectionFilter=NOOPSegmentationPreselection<3>,
        class TRealType=float>
class  RegressionNonLocalFilter : 
    public itk::ImageToImageFilter<TFeatureImage, TOutputImage>
  {
  public:
    
    /** Standard class typedefs. */
    typedef RegressionNonLocalFilter Self;
    typedef itk::ImageToImageFilter<TFeatureImage, TOutputImage>  Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;
    
  
    /** Standard New method. */
    itkNewMacro(Self);  

    /** Runtime information support. */
    itkTypeMacro(RegressionNonLocalFilter,itk::ImageToImageFilter);
    
    /** Image related typedefs. */
    typedef TFeatureImage                                   FeatureImageType;
    typedef TLabelImage                                     LabelImageType;
    typedef TOutputImage                                    OutputImageType;
    
    typedef typename TFeatureImage::RegionType              RegionType;
    typedef typename TFeatureImage::SizeType                SizeType;
    typedef typename TFeatureImage::IndexType               IndexType;
    typedef typename TFeatureImage::PixelType               PixelType;
    typedef typename TLabelImage::PixelType                 LabelPixelType;
    typedef typename TOutputImage::PixelType                OutputPixelType;
    
    typedef typename Superclass::OutputImageRegionType      OutputImageRegionType;
    typedef typename TOutputImage::Pointer            		  OutputImagePointerType;
    
    typedef typename std::vector<ImagePair<TFeatureImage,TLabelImage> > SegmentationLibraryType;
    
#if (ITK_VERSION_MAJOR==3)  
    typedef typename itk::OrientedImage<unsigned char, TFeatureImage::ImageDimension> RoiImageType;
#else
    typedef typename itk::Image<unsigned char, TFeatureImage::ImageDimension> RoiImageType;
#endif    
    typedef typename RoiImageType::Pointer                RoiImagePointer;

    
    /** Image related typedefs. */
    itkStaticConstMacro(ImageDimension, unsigned int, TFeatureImage::ImageDimension);
    typedef typename std::vector<LabelDistance<LabelPixelType,ImageDimension> > SegmentationLabelDistance;

    
    /** Parameters typedefs */
    typedef typename TRegressor::Pointer          RegressorPointer;
    typedef typename TPreselectionFilter::Pointer PreselectionFilterPointerType;
    
    /** Typedef for boundary conditions. */
    typedef itk::ImageBoundaryCondition<FeatureImageType> *       ImageBoundaryConditionPointerType;
    typedef itk::ImageBoundaryCondition<FeatureImageType> const * ImageBoundaryConditionConstPointerType;
    typedef itk::ConstantBoundaryCondition<FeatureImageType>      DefaultBoundaryConditionType;
    
    typedef itk::ImageBoundaryCondition<LabelImageType> *       LabelImageBoundaryConditionPointerType;
    typedef itk::ImageBoundaryCondition<LabelImageType> const * LabelImageBoundaryConditionConstPointerType;
    typedef itk::ConstantBoundaryCondition<LabelImageType>      DefaultLabelBoundaryConditionType;
    

  /** Neighborhood iterator type. */
    typedef itk::ConstNeighborhoodIterator<TFeatureImage>   NeighborhoodIteratorType;
    
    typedef itk::ConstNeighborhoodIterator<TLabelImage>     LabelNeighborhoodIteratorType;

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
    
    /** RegressionNonLocalFilter need to make sure they request enough of an
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
    

    void SetSegmentationLibrary(SegmentationLibraryType &s)
    {
      m_SegmentationLibrary=s;
    }
    
    SegmentationLibraryType & GetSegmentationLibrary(void)
    {
      return m_SegmentationLibrary;
    }


    /** Set distance functor. */
    itkSetMacro(Regressor, RegressorPointer);
    /** Get distance functor. */
    itkGetMacro(Regressor, RegressorPointer);
    
    /** assign ROI Image */
    itkSetMacro(RoiImage,RoiImagePointer);
    /** get ROI Image */
    itkGetMacro(RoiImage,RoiImagePointer);
    
    /** assign preselection filter */
    itkSetMacro(PreselectionFilter,PreselectionFilterPointerType);
    /** get preselection filter */
    itkGetMacro(PreselectionFilter,PreselectionFilterPointerType);

  
    /** get confidence results */
    itkGetMacro(Confidence,FloatImagePointerType);
    /** set confidence results */
    itkSetMacro(Confidence,FloatImagePointerType);
    
    /** get distance results */
    itkGetMacro(SearchDistance,FloatImagePointerType);
    /** set distance results */
    itkSetMacro(SearchDistance,FloatImagePointerType);

    /** get Grading map */
    itkGetMacro(GradingMap,FloatImagePointerType);
    /** set Grading map */
    itkSetMacro(GradingMap,FloatImagePointerType);
    
    /** assign confidence threshold */
    itkSetMacro(ConfidenceThreshold,TRealType);
    /** get confidence threshold */
    itkGetMacro(ConfidenceThreshold,TRealType);

    /** assign label count */
    itkSetMacro(LabelCount,int);
    /** get label count */
    itkGetMacro(LabelCount,int);
    
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
    
    FloatImagePointerType GetProbabilityMap(int i)
    {
      return m_ProbabilityMaps[i];
    }
  protected:
    RegressionNonLocalFilter();
    ~RegressionNonLocalFilter() {};

    /** prepare before starting the thread */
    void BeforeThreadedGenerateData(void) override;

    /** Multi-thread version GenerateData. */
    void  ThreadedGenerateData (const OutputImageRegionType& outputRegionForThread,
                                itk::ThreadIdType threadId) override;
    /** finish calculations */
    void AfterThreadedGenerateData(void) override;
    
    

    void PrintSelf(std::ostream& os, itk::Indent indent) const;

  private:
    RegressionNonLocalFilter(const Self&); //purposely not implemented
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
    /** Default label boundary condition */
    DefaultLabelBoundaryConditionType 
                            m_DefaultLabelBoundaryCondition;
    /** Label boundary condition */
    LabelImageBoundaryConditionPointerType 
                            m_LabelBoundaryCondition;
    /** Regression functor */
    RegressorPointer        m_Regressor;
    
    /** Number of labels including background*/
    int                     m_LabelCount;
    
    /** Segmentation library */
    SegmentationLibraryType m_SegmentationLibrary;
    /** ROI mask */
    RoiImagePointer         m_RoiImage;
    /** Preselection filter */
    PreselectionFilterPointerType     m_PreselectionFilter;
    /** Voxel-level confidence value */
    FloatImagePointerType   m_Confidence;
    /** Voxel-level mean weighted search distance */
    FloatImagePointerType   m_SearchDistance;
    /** Grading map*/
    FloatImagePointerType   m_GradingMap;
    /** Probability maps */
    std::vector<FloatImagePointerType> m_ProbabilityMaps;
    /** Use threshold on confidence map to determine where algorithm have to work, needed for interative approach*/
    TRealType               m_ConfidenceThreshold;
    /** Number of voxels which had confidence below threshold in last iteration */
    TRealType               m_NonConfidentVoxels;
    /** Number of voxels which had confidence below threshold in last iteration (in one thread)*/
    std::vector<TRealType>  m_ThreadNonConfidentVoxels;
    
    /** preloaded results, to initialize segmentation */
    OutputImagePointerType  m_PreloadedResults;
  }; // end of class

} //minc


#include "itkRegressionNonLocalFilter.txx"

#endif //__mincRegressionNonLocalFilter_h

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;show-tabs on;tab-width 2;hl c++
