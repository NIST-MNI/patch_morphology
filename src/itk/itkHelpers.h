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
#ifndef __itkHelpers_h__
#define __itkHelpers_h__

#include <itkMetaDataObject.h>
#include <itkImage.h>
#include <itkExceptionObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>


#include <stdlib.h>
#include <time.h>

//helper functions to simplify common operations
#define REPORT_ERROR(MSG) throw itk::ExceptionObject(__FILE__,__LINE__,MSG)

namespace itk
{
  //! allocate volume of the same dimension,spacing and origin
  template<class T,class S> void allocate_same(typename T::Pointer &image,const typename S::Pointer &sample)
  {
    image->SetLargestPossibleRegion(sample->GetLargestPossibleRegion());
    image->SetBufferedRegion(sample->GetLargestPossibleRegion());
    image->SetRequestedRegion(sample->GetLargestPossibleRegion());
    image->SetSpacing( sample->GetSpacing() );
    image->SetOrigin ( sample->GetOrigin() );
    image->SetDirection(sample->GetDirection());
    image->Allocate();
  }
  
  //! allocate volume of the same dimension,spacing and origin
  template<class T,class S> void allocate_same(typename T::Pointer &image,const typename S::ConstPointer &sample)
  {
    image->SetLargestPossibleRegion(sample->GetLargestPossibleRegion());
    image->SetBufferedRegion(sample->GetLargestPossibleRegion());
    image->SetRequestedRegion(sample->GetLargestPossibleRegion());
    image->SetSpacing( sample->GetSpacing() );
    image->SetOrigin ( sample->GetOrigin() );
    image->SetDirection(sample->GetDirection());
    image->Allocate();
  }
  
  //! allocate volume of the same dimension,spacing and origin
  template<class T,class S> typename T::Pointer allocate_same(const typename S::Pointer &sample)
  {
    typename T::Pointer image=T::New();
    image->SetLargestPossibleRegion(sample->GetLargestPossibleRegion());
    image->SetBufferedRegion(sample->GetLargestPossibleRegion());
    image->SetRequestedRegion(sample->GetLargestPossibleRegion());
    image->SetSpacing( sample->GetSpacing() );
    image->SetOrigin ( sample->GetOrigin() );
    image->SetDirection(sample->GetDirection());
    image->Allocate();
    return image;
  }

  static void copy_metadata(itk::Object* dst,itk::Object* src)
  {
    dst->SetMetaDataDictionary(src->GetMetaDataDictionary());
  }

  template <class T> typename T::Pointer load_volume(const char *file)
  {
     
    typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();
    
    reader->SetFileName(file);
    reader->Update();
    
    return reader->GetOutput();
  }

  template <class T> void save_volume(const char *file,typename T::Pointer img)
  {
    typename itk::ImageFileWriter< T >::Pointer writer = itk::ImageFileWriter<T>::New();
    writer->SetFileName(file);
    writer->SetInput( img );
    
    if( getenv("MINC_COMPRESS") != NULL)
      writer->SetUseCompression( true );
    
    writer->Update();
  } 
  
  
  // MINC file specific helpers
  static void append_minc_history(itk::Object* dst,const std::string& history)
  {
    std::string old_history;
    itk::ExposeMetaData<std::string>( dst->GetMetaDataDictionary() , "history",old_history);
    old_history+=history;
    itk::EncapsulateMetaData( dst->GetMetaDataDictionary(),"history",old_history);
  }
  
  static void set_minc_storage_type(itk::Object* image,std::string datatype)
  {
    itk::EncapsulateMetaData<std::string>(image->GetMetaDataDictionary(),"storage_data_type",datatype);
  }
  
  static void copy_minc_dimorder(itk::Object* dst,itk::Object* src)
  {
    std::string dimorder;
    if(itk::ExposeMetaData(src->GetMetaDataDictionary(),"dimension_order",dimorder))
    {
      itk::EncapsulateMetaData(dst->GetMetaDataDictionary(),"dimension_order",dimorder);
    }
  }

  //! default label voxel type
  typedef unsigned char mask_voxel;
  
  //! default minc file voxel type
  typedef float voxel_type;
  
  //! default minc volume dimension
  const int volume_dimensions = 3;

  typedef itk::Image < voxel_type, volume_dimensions > image3d;
  typedef itk::Image < mask_voxel, volume_dimensions > mask3d;

  
  static std::string minc_time_stamp(int argc,char **argv)
  {
    std::string timestamp;
    
    char cur_time[200];
    time_t t;
    struct tm *tmp;

    t = time(NULL);
    tmp = localtime(&t);
    
    strftime(cur_time, sizeof(cur_time), "%a %b %d %T %Y>>>", tmp);
    /* Get the time, overwriting newline */
    timestamp=cur_time;
    
    /* Copy the program name and arguments */
    for (int i=0; i<argc; i++) {
      timestamp+=argv[i];
      timestamp+=" ";
    }
    timestamp+="\n";
    
    return timestamp;
  }
}

#endif //__itkHelpers_h__
// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;show-tabs on;tab-width 2;hl c++


