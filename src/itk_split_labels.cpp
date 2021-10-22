/* ----------------------------- MNI Header -----------------------------------
@NAME       :  itk_split_labels
@DESCRIPTION:  split discrete labels into individual volumes with antialiasing
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
#include <math.h>

#include <itkVector.h>
#include <itkResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkVectorImage.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageConstIterator.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkAntiAliasBinaryImageFilter.h>
#include <itkExpNegativeImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkMultiplyImageFilter.h>

#include <itkNearestNeighborInterpolateImageFunction.h>

#include <unistd.h>
#include <getopt.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>
#include <itkImageIOBase.h>

#include <vnl/vnl_cross.h>

#include "itkHelpers.h"

typedef itk::ImageBase<3>           Image3DBase;
typedef itk::Image<float,3>         Float3DImage;
typedef itk::Image<int,3>           Int3DImage;
typedef itk::Image<short,3>         Short3DImage;
typedef itk::Image<unsigned char,3> Byte3DImage;

typedef itk::ImageIOBase          IOBase;
typedef itk::SmartPointer<IOBase> IOBasePointer;

typedef itk::BSplineInterpolateImageFunction< Float3DImage, double, double >  InterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction< Float3DImage, double >    NNInterpolatorType;

typedef itk::ResampleImageFilter<Float3DImage, Float3DImage> FloatFilterType;


using namespace  std;

void show_usage (const char * prog)
{
  std::cerr
    << "Usage: "<<prog<<" <input> <output_XXXX.mnc>  where XXXX is printf format string" << std::endl
    << "--clobber overwrite files"    << std::endl
    << "--order <n> spline order, default 2 "<<std::endl
    << "--uniformize <step> - will make a volume with uniform step size and no direction cosines" << std::endl
    << "--unistep <step> - will make a volume with the same step size and preserve direction cosines" << std::endl
    << "--byte  - store image in byte  voxels minc file"<<std::endl
    << "--short - store image in short voxels minc file"<<std::endl
    << "--float - store image in float voxels minc file"<<std::endl
    << "--relabel map.txt apply relabeling map"<<std::endl
    << "--lut map.txt  apply relabeling map"<<std::endl
    << "--lut-string \"a b;c d;....\" apply lut string in command line"<<std::endl
    << "--aa apply binary anti-aliasing filter "<<std::endl
    << "--antialias apply binary anti-aliasing filter "<<std::endl
    << "--blur <fwhm> apply blurring filter "<<std::endl
    << "--expit <beta> apply expit (inverse-logit) 1/(1+exp(-x*beta)) function to output "<<std::endl
    << "--normalize normalize probabiblites for expit output"<<std::endl
    << "--missing <max> output missing labels up to max id"<<std::endl
    << "--layers <n> specify maximum number of layers (distance) from the boundary, default 4"<<std::endl;
}

template<class T,class I> void generate_uniform_sampling(T* flt, const I* img,double step)
{
  //obtain physical coordinats of all image corners
  typename I::RegionType r=img->GetLargestPossibleRegion();
  std::vector<double> corner[3];
  for(int i=0;i<8;i++)
  {
    typename I::IndexType idx;
    typename I::PointType c;
    idx[0]=r.GetIndex()[0]+r.GetSize()[0]*(i%2);
    idx[1]=r.GetIndex()[1]+r.GetSize()[1]*((i/2)%2);
    idx[2]=r.GetIndex()[2]+r.GetSize()[2]*((i/4)%2);
    img->TransformIndexToPhysicalPoint(idx,c);
    for(int j=0;j<3;j++)
      corner[j].push_back(c[j]);
  }
  typename I::IndexType start;
  typename T::SizeType size;
  typename T::OriginPointType org;
  typename I::SpacingType spc;
  spc.Fill(step);
  for(int j=0;j<3;j++)
  {
    std::sort(corner[j].begin(),corner[j].end());
    size[j]=ceil((corner[j][7]-corner[j][0])/step);
    org[j]=corner[j][0];
  }
  Float3DImage::DirectionType identity;
  identity.SetIdentity();
  flt->SetOutputDirection(identity);
  start.Fill(0);
  flt->SetOutputStartIndex(start);
  flt->SetSize(size);
  flt->SetOutputOrigin(org);
  flt->SetOutputSpacing(spc);
}

template<class T,class I> void generate_unistep_sampling(T* flt, const I* img,double step)
{
  //obtain physical coordinats of all image corners
  typename I::RegionType r=img->GetLargestPossibleRegion();

  typename I::IndexType start;
  typename T::SizeType size;
  typename T::OriginPointType org=img->GetOrigin();
  typename I::SpacingType spc;
  spc.Fill(step);

  for(int j=0;j<3;j++)
  {
    org[j]=org[j]-img->GetSpacing()[j]/2.0+step/2.0;
    size[j]=::ceil((r.GetSize()[j]+1)*img->GetSpacing()[j]/step);
  }

  start.Fill(0);

  flt->SetOutputDirection(img->GetDirection());

  flt->SetOutputStartIndex(start);
  flt->SetSize(size);
  flt->SetOutputOrigin(org);
  flt->SetOutputSpacing(spc);
}


template<class Image,class ImageOut,class Interpolator>
void resample_label_image (
   IOBase* base,
   const std::string& output_f,
   double uniformize,
   double unistep,
   const char* history,
   bool store_float,
   bool store_short,
   bool store_byte,
   Interpolator* interpolator,
   const std::map<int,int>& label_map=std::map<int,int>(),
   double blur_fwhm=0.0,
   bool   baa_smooth=false,
   double expit_beta=0.0,
   bool   normalize_p=false,
   int    output_missing=0,
   int    layers=4
                          )
{
  typedef typename itk::ResampleImageFilter<ImageOut, ImageOut> ResampleFilterType;
  typedef typename itk::ImageFileReader<Image >                 ImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >              ImageWriterType;
  
  typedef itk::ImageRegionConstIterator<Image>                  ConstInputImageIteratorType;
  typedef itk::ImageRegionConstIterator<ImageOut>               ConstTmpImageIteratorType;
  
  typedef itk::BinaryThresholdImageFilter<Image,ImageOut>       ThresholdFilterType;
  typedef itk::SmoothingRecursiveGaussianImageFilter<ImageOut,ImageOut>  BlurFilterType;
  typedef itk::AntiAliasBinaryImageFilter<ImageOut,ImageOut>    BAAFilterType;
  typedef itk::ExpNegativeImageFilter< ImageOut, ImageOut > InvExpFilterType;  
  typedef itk::AddImageFilter< ImageOut, ImageOut > AddFilterType;
  typedef itk::DivideImageFilter< ImageOut,ImageOut, ImageOut > DivideFilterType;
  typedef itk::ImageDuplicator< ImageOut > ImageDuplicatorType;
  typedef itk::MultiplyImageFilter< ImageOut,ImageOut, ImageOut > MultiplyFilterType;

  typedef typename Image::PixelType InputPixelType;
  typename ImageReaderType::Pointer reader = ImageReaderType::New();

  double blur_sigma=blur_fwhm/(2.0*sqrt(2*log(2.0)));
  bool do_resample=false;

  //initializing the reader
  reader->SetImageIO(base);
  reader->SetFileName(base->GetFileName());
  reader->Update();

  typename Image::Pointer in=reader->GetOutput();

  if(!label_map.empty())
  {
    for(itk::ImageRegionIterator<Image> it(in, in->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
      std::map<int,int>::const_iterator pos=label_map.find(static_cast<int>(it.Get()));
      if(pos==label_map.end())
        it.Set(0); //set to BG
      else
        it.Set(static_cast<InputPixelType>((*pos).second));
    }
  }

  typename ResampleFilterType::Pointer  filter           = ResampleFilterType::New();
  typename ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
  typename BlurFilterType::Pointer      blur_filter      = BlurFilterType::New();
  typename BAAFilterType::Pointer       baa_filter       = BAAFilterType::New();

  typename InvExpFilterType::Pointer   inv_exp_filter    = InvExpFilterType::New();
  typename AddFilterType::Pointer      add_filter        = AddFilterType::New();
  typename AddFilterType::Pointer      sum_filter        = AddFilterType::New();
  typename DivideFilterType::Pointer   divide_filter     = DivideFilterType::New();
  typename MultiplyFilterType::Pointer multiply_filter   = MultiplyFilterType::New();
  

  typename ImageDuplicatorType::Pointer duplicator       = ImageDuplicatorType::New();

  //creating the interpolator
  filter->SetInterpolator( interpolator );
  filter->SetDefaultPixelValue( 0 );

  //this is for processing using batch system
  filter->SetNumberOfThreads(1);

  typename Image::Pointer like=nullptr;

  if(uniformize!=0.0)
  {
    generate_uniform_sampling<ResampleFilterType,Image>(filter,in,uniformize);
    do_resample=true;
  } else if(unistep!=0.0) {
    generate_unistep_sampling<ResampleFilterType,Image>(filter,in,unistep);
    do_resample=true;
  } else {
    //we are using original sampling
    filter->SetOutputParametersFromImage(in);
    filter->SetOutputDirection(in->GetDirection());
  }
    
  typename ImageOut::RegionType region;
  region.SetSize (filter->GetSize());
  region.SetIndex(filter->GetOutputStartIndex());

  //split the input image
  std::set<InputPixelType> sval;
  for(ConstInputImageIteratorType it(in, in->GetBufferedRegion()); !it.IsAtEnd(); ++it)
  {
    InputPixelType val = it.Get();

    if(isfinite(val))
      sval.insert(val);
  }

  std::vector<typename ImageOut::Pointer> output_images;
  
  //TODO: process labels in parallel here?
  for(typename set<InputPixelType>::iterator label_it = sval.begin(); label_it != sval.end(); ++label_it)
  {
    threshold_filter->SetLowerThreshold(*label_it);
    threshold_filter->SetUpperThreshold(*label_it);
    threshold_filter->SetInsideValue(1.0);
    threshold_filter->SetOutsideValue(0.0);
    threshold_filter->SetInput(in);

    threshold_filter->Update();
    typename ImageOut::Pointer img=threshold_filter->GetOutput();
    
    if(baa_smooth) {
      //baa_filter->SetMaximumRMSChange(baa_smooth);
      baa_filter->SetNumberOfLayers(layers);
      baa_filter->SetInput(img);
      baa_filter->Update();
      img=baa_filter->GetOutput();
    }
    
    if(blur_fwhm>0.0)
    {
      //blur_filter->SetUseImageSpacing(true);
      blur_filter->SetSigma(blur_sigma);
      blur_filter->SetInput(img);
      blur_filter->Update();
      img=blur_filter->GetOutput(); 
    }

    if(do_resample)
    {
      filter->SetInput(img);
      filter->Update();
      img=filter->GetOutput();
    }

    if( expit_beta > 0.0 )
    {
      multiply_filter->SetInput1(img);
      multiply_filter->SetConstant2(expit_beta);
      
      inv_exp_filter->SetInput(multiply_filter->GetOutput());
      
      add_filter->SetInput1(inv_exp_filter->GetOutput());
      add_filter->SetConstant2(1.0);

      divide_filter->SetConstant1(1.0);
      divide_filter->SetInput2(add_filter->GetOutput());

      divide_filter->Update();
      
      img=divide_filter->GetOutput();
    }
    
    duplicator->SetInputImage( img );
    duplicator->Update();
    img=duplicator->GetOutput();

    output_images.push_back(img);


    if( normalize_p && expit_beta > 0.0 )
    {
      if(label_it == sval.begin())
      {
        sum_filter->SetInput1(img);
      }
      else
      {
        sum_filter->SetInput2(img);
        sum_filter->Update();
        sum_filter->SetInput1(sum_filter->GetOutput());
      }
    }
  }

  int i=0;
  
  for(typename set<InputPixelType>::iterator label_it = sval.begin(); label_it != sval.end(); ++label_it, ++i)
  {
    typename ImageOut::Pointer img=output_images[i];

    if( normalize_p && expit_beta > 0.0 )
    {
      sum_filter->SetConstant2(0.0);
      sum_filter->Update();
      divide_filter->SetInput1(img);
      divide_filter->SetInput2(sum_filter->GetOutput());
      divide_filter->Update();
      img=divide_filter->GetOutput();
    }

    itk::copy_metadata(img, in);
    itk::append_minc_history(img, history);
    
    //generic file writer
    typename ImageWriterType::Pointer writer = ImageWriterType::New();
    char label_file[1024];
    sprintf(label_file, output_f.c_str(), (int)(*label_it));
    
    writer->SetFileName(label_file);
    if(store_float)
    {
      itk::set_minc_storage_type(img,typeid(float).name());
    } else if(store_short) {
      itk::set_minc_storage_type(img,typeid(unsigned short).name());
    } else if(store_byte) {
      itk::set_minc_storage_type(img,typeid(unsigned char).name());
    }
    writer->SetInput( img );

    if( getenv("MINC_COMPRESS") != NULL)
      writer->SetUseCompression( true );

    writer->Update();
    std::cout<<(int)*label_it<<","<<label_file<<std::endl;
  }
  
  if(output_missing>0) // now let's output missing empty labels
  {
    typename ImageOut::Pointer img=output_images[0];
    img->FillBuffer(0.0);
    
    for(int i=0;i<output_missing;i++)
    {
      if(sval.find(i)==sval.end())
      {
        typename ImageWriterType::Pointer writer = ImageWriterType::New();
        char label_file[1024];
        sprintf(label_file, output_f.c_str(), i);
    
        writer->SetFileName(label_file);
        if(store_float)
        {
          itk::set_minc_storage_type(img,typeid(float).name());
        } else if(store_short) {
          itk::set_minc_storage_type(img,typeid(unsigned short).name());
        } else if(store_byte) {
          itk::set_minc_storage_type(img,typeid(unsigned char).name());
        }
        writer->SetInput( img );

        if( getenv("MINC_COMPRESS") != NULL)
          writer->SetUseCompression( true );

        writer->Update();
        std::cout<<i<<","<<label_file<<std::endl;
      }
    }
  }
  
  
}



int main (int argc, char **argv)
{
  int store_float=0;
  int store_short=0;
  int store_byte=0;

  int verbose=0, clobber=0;
  int order=2;
  double   uniformize=0.0;
  double   unistep=0.0;
  int      invert=0;
  int      labels=0;
  std::string       history;
  std::string       map_f,map_str;
  std::map<int,int> label_map;
  int    baa_smooth=0;
  double blur_fwhm=0.0;
  double expit_beta=0.0;
  int    normalize_p=0;
  int    output_missing=0;
  int    layers=4;

  history = itk::minc_time_stamp(argc,argv);
  bool order_was_set=false;

  static struct option long_options[] = {
    {"verbose",    no_argument,       &verbose, 1},
    {"quiet",      no_argument,       &verbose, 0},
    {"clobber",    no_argument,       &clobber, 1},
    {"missing",    required_argument, 0 , 'm'},
    {"order",      required_argument, 0, 'o'},
    {"uniformize", required_argument, 0, 'u'},
    {"unistep",    required_argument, 0, 'U'},
    {"float",      no_argument,       &store_float, 1},
    {"short",      no_argument,       &store_short, 1},
    {"byte",       no_argument,       &store_byte,  1},
    {"relabel",    required_argument,      0,'L'},
    {"lut",        required_argument,      0,'L'},
    {"lut-string", required_argument,      0,'s'},
    {"blur",       required_argument,      0,'b'},
    {"antialias",  no_argument,       &baa_smooth, 1},
    {"aa",         no_argument,       &baa_smooth, 1},
    {"expit",      required_argument, 0, 'E'},
    {"normalize",  no_argument,       &normalize_p, 1},
    {"layers",     required_argument,      0,'l'},
    {0, 0, 0, 0}
    };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "vqcl:t:o:u:L:l:s:a:U:b:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;

    switch (c)
    {
    case 0:
        break;
    case 'v':
        cout << "Version: 0.8" << endl;
        return 0;
    case 'o':
      order=atoi(optarg);
      order_was_set=true;
      break;
    case 'u':
      uniformize=atof(optarg);break;
    case 'U':
      unistep=atof(optarg);break;
    case 'b':
      blur_fwhm=atof(optarg);break;
    case 'E':
      expit_beta=atof(optarg);break;
    case 'l':
      layers=atoi(optarg);break;
    case 'L':
      map_f=optarg;break;
    case 's':
      map_str=optarg;break;
    case 'm':
      output_missing=atoi(optarg);break;
    case '?':
        /* getopt_long already printed an error message. */
    default:
        show_usage (argv[0]);
        return 1;
    }
  }
  
  if ((argc - optind) < 2) {
      show_usage(argv[0]);
      return 1;
  }

  std::string input_f=argv[optind];
  std::string output_f=argv[optind+1];

//   if (!clobber && !access (output_f.c_str (), F_OK))
//   {
//     std::cerr << output_f.c_str () << " Exists!" << std::endl;
//     return 1;
//   }

  if(!map_f.empty())
  {
    if(verbose)
      std::cout<<"Going to relabel input according to:"<<map_f.c_str()<<std::endl;
    std::ifstream in_map(map_f.c_str());

    while(!in_map.eof() && in_map.good())
    {
      int l1,l2;
      in_map>>l1>>l2;
      std::map<int,int>::iterator pos=label_map.find(l1);
      if(pos==label_map.end())
        label_map[l1]=l2;
      else if(label_map[l1]!=l2 and verbose)
      {
        std::cerr<<"Warning: label "<<l1<<" mapped more then once!"<<std::endl;
        std::cerr<<"Previous map:"<<(*pos).second<<" New map:"<<l2<<std::endl;
      }
    }
  }

  if(!map_str.empty())
  {
    if(verbose)
      std::cout<<"Going to relabel input according to:"<<map_str.c_str()<<std::endl;

    const char* delim1=";";
    char *saveptr1;
    char *_str=strdup(map_str.c_str());

    for(char *tok1=strtok_r(_str,delim1,&saveptr1);tok1;tok1=strtok_r(NULL,delim1,&saveptr1))
    {
      int l1,l2;
      if(sscanf(tok1,"%d %d",&l1,&l2)!=2)
      {
        std::cerr<<"Unrecognized lut string:"<<tok1<<std::endl;
      } else {
        std::map<int,int>::iterator pos=label_map.find(l1);
        if(pos==label_map.end())
          label_map[l1]=l2;
        else if(label_map[l1]!=l2 and verbose)
        {
          std::cerr<<"Warning: label "<<l1<<" mapped more then once!"<<std::endl;
          std::cerr<<"Previous map:"<<(*pos).second<<" New map:"<<l2<<std::endl;
        }
      }
    }
    free(_str);
  }

  try
  {
    //try to figure out what we have got
    IOBasePointer io = itk::ImageIOFactory::CreateImageIO(input_f.c_str(), itk::ImageIOFactory::ReadMode );

    if(!io)
      throw itk::ExceptionObject("Unsupported image file type");

    io->SetFileName(input_f.c_str());
    io->ReadImageInformation();

    size_t nd = io->GetNumberOfDimensions();
    size_t nc = io->GetNumberOfComponents();

    if(verbose)
    {
      if(uniformize!=0.0)
        std::cout<<"Making uniform sampling, step size="<<uniformize<<std::endl;
      else if(unistep)
        std::cout<<"Making same step size, new step size="<<unistep<<std::endl;

      if(normalize_p)
        std::cout<<"Performing probability normalization"<<std::endl;
    }

    if( nc==1 && nd==3 ) //3D image, simple case
    {
          InterpolatorType::Pointer interpolator = InterpolatorType::New();
          interpolator->SetSplineOrder(order);

          resample_label_image<Int3DImage,Float3DImage,InterpolatorType>(
            io,output_f,uniformize,unistep,
            history.c_str(),
            store_float,store_short,store_byte,
            interpolator,
            label_map,blur_fwhm,baa_smooth,expit_beta,normalize_p,output_missing,layers
          );

    } else {
      throw itk::ExceptionObject("This number of dimensions is not supported currently");
    }
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

// kate: space-indent on; indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;show-tabs on
