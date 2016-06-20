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
#ifndef __mincNonLocalPatchesFilter_h
#define __mincNonLocalPatchesFilter_h


#include "itkImageToImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhood.h"
#include "itkConstSliceIterator.h"
#include "itkImageBoundaryCondition.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImageRegionIterator.h"

#include "itkPreselectionFilter.h"


// FOR Debug purposes only

//#include "itkMincHelpers.h"

namespace itk {

/** \class NonLocalPatchesFilter 
 * \brief Base class for the non-local patch-based filters
 *
 *
 * based on itk::MorphologyImageFilter filter
 * 
 * Subclasses of this class can define their own operations by simply
 * providing their own Evaluate() protected member function.
 *
 */

template<class TInputImage, 
         class TOutputImage, 
         class TSearch,
         class TPatch,
         class TPreselectionFilter=NOOPPreselection<3> >
class NonLocalPatchesFilter : 
    public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard Self typedef */
  typedef NonLocalPatchesFilter                              Self;
  typedef itk::ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef itk::SmartPointer<Self>                            Pointer;
  typedef itk::SmartPointer<const Self>                      ConstPointer;
  
  /** Runtime information support. */
  itkTypeMacro(NonLocalPatchesFilter, itk::ImageToImageFilter);

  /** Standard New method. */
  itkNewMacro(Self);  
  
  
  /** Image related typedefs. */
  typedef TInputImage                                   InputImageType;
  typedef TOutputImage                                  OutputImageType;
  typedef typename TInputImage::RegionType              RegionType;
  typedef typename TInputImage::SizeType                SizeType;
  typedef typename TInputImage::IndexType               IndexType;
  typedef typename TInputImage::PixelType               PixelType;
  typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;
  typedef typename TOutputImage::Pointer            		OutputImagePointer;
#if (ITK_VERSION_MAJOR==3)
  typedef typename itk::OrientedImage<unsigned char, TInputImage::ImageDimension> RoiImageType;
#else
  typedef typename itk::Image<unsigned char, TInputImage::ImageDimension> RoiImageType;
#endif

  typedef typename RoiImageType::Pointer                RoiImagePointer;

  
  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Typedef for boundary conditions. */
  typedef itk::ImageBoundaryCondition<InputImageType> *       ImageBoundaryConditionPointerType;
  typedef itk::ImageBoundaryCondition<InputImageType> const * ImageBoundaryConditionConstPointerType;
  typedef itk::ConstantBoundaryCondition<InputImageType>      DefaultBoundaryConditionType;
  

/** Neighborhood iterator type. */
  typedef itk::ConstNeighborhoodIterator<TInputImage>  	NeighborhoodIteratorType;

  /** Kernel typedef. */
  typedef TSearch SearchKernelType;
  typedef TPatch  PatchKernelType;
  
  typedef typename TPreselectionFilter::Pointer PreselectionFilterPointerType;
  
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
  
  /** NonLocalPatchesFilter need to make sure they request enough of an
   * input image to account for the structuring element size.  The input
   * requested region is expanded by the radius of the Search structuring element.
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
  
  /** Get the current boundary condition. */
  itkGetConstMacro(BoundaryCondition, ImageBoundaryConditionPointerType);

	/** assign ROI Image */
	itkSetMacro(RoiImage,RoiImagePointer);
	/** get ROI Image */
	itkGetMacro(RoiImage,RoiImagePointer);
	
	/** assign preselection filter */
	itkSetMacro(PreselectionFilter,PreselectionFilterPointerType);
	/** get preselection filter */
	itkGetMacro(PreselectionFilter,PreselectionFilterPointerType);
	
	
	/** for debugging */
	itkGetMacro(Weights, OutputImagePointer);
	itkGetMacro(SmallestWeights, OutputImagePointer);

	
  
protected:
  NonLocalPatchesFilter();
  ~NonLocalPatchesFilter() {};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;
	
	/** prepare before starting the thread */
	void BeforeThreadedGenerateData(void);

#if (ITK_VERSION_MAJOR==3)
  /** Multi-thread version GenerateData. */
  void  ThreadedGenerateData (const OutputImageRegionType& ,
                              int threadId);
#else
  /** Multi-thread version GenerateData. */
  void  ThreadedGenerateData (const OutputImageRegionType & ,
                              itk::ThreadIdType threadId);
#endif
  
  /** Evaluate image neighborhood with kernel to find the new value 
   * for the center pixel value. */
  virtual PixelType Evaluate(const NeighborhoodIteratorType &searchIt,
                              const NeighborhoodIteratorType &patchIt1,
                              NeighborhoodIteratorType &patchIt2,
                              const SearchKernelIteratorType searchKernelBegin,
                              const SearchKernelIteratorType searchKernelEnd,
                              const PatchKernelIteratorType patchKernelBegin,
                              const PatchKernelIteratorType patchKernelEnd,
                              PreselectionFilterPointerType flt
                            )
  {
    itkExceptionMacro( << "Not Implemented" );
  }

private:
  NonLocalPatchesFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

	
protected:	
  /** kernel or structuring element to use. */
  SearchKernelType m_SearchKernel;
	
  /** kernel or structuring element to use. */
  PatchKernelType m_PatchKernel;

  /** Pointer to a persistent boundary condition object used
   * for the image iterator. */
  ImageBoundaryConditionPointerType m_BoundaryCondition;

  /** Default boundary condition */
  DefaultBoundaryConditionType   m_DefaultBoundaryCondition;
	
	OutputImagePointer             m_Weights;
	OutputImagePointer             m_SmallestWeights;
	RoiImagePointer                m_RoiImage;
	PreselectionFilterPointerType  m_PreselectionFilter;
	
}; // end of class

} // end namespace itk
  
//#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNonLocalPatchesFilter.txx"
//#endif

#endif

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2
