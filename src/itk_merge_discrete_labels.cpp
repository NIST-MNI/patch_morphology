/* ----------------------------- MNI Header -----------------------------------
@NAME       :  itk_resample
@DESCRIPTION:  an example of using spline itnerpolation with MINC xfm transform
@COPYRIGHT  :
              Copyright 2011 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif 

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <valarray>
#include <math.h>
#include <limits>
#include <unistd.h>
#include <algorithm>
#include <stdlib.h>

#include <itkIdentityTransform.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageConstIterator.h>

#include <unistd.h>
#include <getopt.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>
#include <itkImageIOBase.h>

#include <vnl/vnl_cross.h>

#include "itkHelpers.h"

typedef itk::Image<int,3>           Int3DImage;
typedef itk::Image<short,3>         Short3DImage;
typedef itk::Image<unsigned char,3> Byte3DImage;

using namespace  std;

void show_usage (const char * prog)
{
  std::cerr 
    << "Usage: "<<prog<<" input_1.mnc [input_2 ....] <output.mnc> " << std::endl
    << "--csv <input> read input from csv file" << std::endl
    << "--clobber overwrite files"    << std::endl
    << "--byte  - store image in byte  voxels minc file"<<std::endl
    << "--short - store image in short voxels minc file"<<std::endl;
 
}


template<class ImageOut>
void assemble_image (const std::vector<std::string>& input_f,
                     const std::string & output_f,
                     const std::string & history,
                     bool store_short,
                     bool store_byte)
{
  typedef typename itk::ImageFileReader<Int3DImage >  ImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >    ImageWriterType;
  
  typedef itk::ImageRegionConstIterator<Int3DImage>   ConstInputImageIteratorType;
  typedef itk::ImageRegionIterator<Int3DImage>        InputImageIteratorType;
  typedef itk::ImageRegionIterator<ImageOut>          ImageOutIteratorType;

  
  typename ImageReaderType::Pointer reader = ImageReaderType::New();

  typename ImageOut::Pointer LabelImage= ImageOut::New();
  
  std::vector<std::string>::const_iterator file_it;
  
  for(file_it=input_f.begin() ; file_it != input_f.end(); ++file_it)
  {
    reader->SetFileName((*file_it).c_str());
    reader->Update();

    typename Int3DImage::Pointer in=reader->GetOutput();
    if(file_it==input_f.begin())
    {
      typename Int3DImage::RegionType region=in->GetLargestPossibleRegion();
      
      LabelImage->SetOrigin(in->GetOrigin());
      LabelImage->SetSpacing(in->GetSpacing());
      LabelImage->SetDirection(in->GetDirection());

      LabelImage->SetLargestPossibleRegion(region);
      LabelImage->SetBufferedRegion(region);
      LabelImage->SetRequestedRegion(region);
      LabelImage->Allocate();
      LabelImage->FillBuffer(0);

      itk::copy_metadata(LabelImage,in);
    }
    
    ConstInputImageIteratorType it_in(in, reader->GetOutput()->GetBufferedRegion());
    ImageOutIteratorType        it_out(LabelImage, LabelImage->GetBufferedRegion());
    
    for(; !it_in.IsAtEnd(); ++it_in,++it_out)
    {
      typename Int3DImage::PixelType val =      it_in.Get();
      
      if(val!=0)
      {
        it_out.Set(val);
      }
    }
  }
  
  itk::append_minc_history(LabelImage,history);
  
  //generic file writer
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(output_f.c_str());
  
  if(store_short) {
    itk::set_minc_storage_type(LabelImage,typeid(unsigned short).name());
  } else if(store_byte) {
    itk::set_minc_storage_type(LabelImage,typeid(unsigned char).name());
  }

  writer->SetInput( LabelImage );
  
  if( getenv("MINC_COMPRESS") != NULL)
    writer->SetUseCompression( true );
  
  writer->Update();
}


int main (int argc, char **argv)
{
  int store_short=0;
  int store_byte=0;
  
  int verbose=0, clobber=0;
  std::string history;
  std::string output_f;
  std::string csv_f;
  
  history = itk::minc_time_stamp(argc,argv);
   
  
  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, 1},
    {"quiet",   no_argument,       &verbose, 0},
    {"clobber", no_argument,       &clobber, 1},
    {"short",   no_argument, &store_short, 1},
    {"byte",    no_argument, &store_byte,  1},
    {"csv ", required_argument,      0,'c'},
    {0, 0, 0, 0}
    };
  
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "vc:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1) break;

      switch (c)
      {
      case 0:
        break;
      case 'v':
        cout << "Version: 1.0" << endl;
        return 0;
      case 'c':
        csv_f=optarg;break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage (argv[0]);
        return 1;
      }
    }

  if ( ((argc - optind) < 1 && csv_f.empty()) || ((argc - optind) < 3 && csv_f.empty()) ) {
    show_usage(argv[0]);
    return 1;
  }
  
  output_f=argv[argc-1];
  argc--;
  
  if (!clobber && !access (output_f.c_str (), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }

  std::vector<std::string> input_f;

  if(!csv_f.empty())
  {
    if(verbose)
      std::cout<<"Going to read information from :"<<csv_f.c_str()<<std::endl;
    
    std::ifstream in_csv(csv_f.c_str());
    
    while(!in_csv.eof() && in_csv.good())
    {
      int l1;
      char tmp[1024],fname[1024];
      in_csv.getline(tmp,sizeof(fname));
      if(!strlen(tmp)) break;
      
      if( sscanf(tmp,"%s",fname)!=1 )
      {
        std::cerr<<"Unexpected line format:\""<<tmp<<"\""<<std::endl;
        return 1;
      }

      input_f.push_back(fname);
    }
  } else {
    for(int i=optind;i<argc;i++)
    {
      input_f.push_back(argv[i]);
    }
  }

  try
  {
    if(store_byte)
      assemble_image<Byte3DImage> (input_f, output_f, history, store_short, store_byte);
    else if(store_short)
      assemble_image<Short3DImage>(input_f, output_f, history, store_short, store_byte);
    else
      assemble_image<Int3DImage>  (input_f, output_f, history, store_short, store_byte);
    return 0;
  } 
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 2;
  }
  return 0;
};
