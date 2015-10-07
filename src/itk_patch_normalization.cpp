/* ----------------------------- MNI Header -----------------------------------
@NAME       :  itk_patch_normalizetion
@DESCRIPTION:  patch-based intensity normalization 
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

#include "itkHelpers.h"

#include <itkBinaryBallStructuringElement.h>
#include <itkFlatStructuringElement.h>
#include <itkMedianImageFilter.h>
#include <itkExpNegativeImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkStatisticsImageFilter.h>

#include "itkNormalizationNonLocalFilter.h"
#include "itkSegmentationPreselectionFilter.h"

#include <iostream>
#include <getopt.h>
#include <libgen.h>
#include <unistd.h>

typedef itk::Image<float,3> FeatureImageType;
typedef itk::Image<unsigned char,3> LabelImageType;

typedef itk::FlatStructuringElement< 3  >  StructuringElementType;

typedef itk::MedianNormalization<float,float >   NormalizationType;

typedef itk::NormalizationNonLocalFilter<
                            FeatureImageType,
                            FeatureImageType,
                            StructuringElementType,
                            StructuringElementType,
                            NormalizationType,
                            itk::L2PatchDistance<FeatureImageType,StructuringElementType> >
                            NormalizationFilterNoPsType;

typedef itk::MedianImageFilter<FeatureImageType, FeatureImageType >  MedianFilterType;
typedef itk::MultiplyImageFilter<FeatureImageType, FeatureImageType, FeatureImageType >  MultiplyFilterType;
typedef itk::ExpNegativeImageFilter<FeatureImageType, FeatureImageType > NegExpFilterType;
typedef itk::AddImageFilter<FeatureImageType, FeatureImageType, FeatureImageType > AddFilterType;

using namespace itk;
using namespace std;

//! vector of strings 
typedef std::vector<std::string> strings;

//! string table (vector of vectors of strings)
typedef std::vector<strings> string_table;

std::string _dirname(const char *file)
{
  char* tmp=new char[strlen(file)+1];
  strcpy(tmp,file);
  std::string r=::dirname(tmp);
  delete[] tmp;
  return r;
}

std::string _basename(const char *file)
{
  char* tmp=new char[strlen(file)+1];
  strcpy(tmp,file);
  std::string r=::basename(tmp);
  delete[] tmp;
  return r;
}

/** read text table, appending currect directory name to the entries up to limit for each line*/
void read_table_n(const char *fn,string_table& tbl,int limit)
{
  tbl.clear();
  //read input file line by line, ignoring empty lines and lines starting with #
  std::ifstream in(fn);
  if(!in.good())
    REPORT_ERROR("Can't open file for reading!");

  std::string dir=_dirname(fn);
  dir+="/";
  while(in.good() && !in.eof())
  {
    char tmp[32000];
    in.getline(tmp,sizeof(tmp));
    if(!strlen(tmp)) continue;
    if(tmp[0]=='#') continue;
    char *str,*saveptr;
    strings _s;
    //std::cout<<tmp<<std::endl;
    int k=0;
    for(str=tmp;;str=NULL)
    {
      const char *token = strtok_r(str, ",", &saveptr);
      if(!token || !strlen(token)) break;
      if(k<limit)
        _s.push_back(dir+token);
      else
        _s.push_back(token);
      k++;
    }
    if(_s.empty()) continue;
    tbl.push_back(_s);
  }
}

bool table_verify(const string_table& tbl)
{
  if(tbl.empty()) 
  {
    std::cerr<<"No data"<<std::endl;
    return false; //?
  }
  int sz=tbl[0].size();
  if(!sz)
  {
    std::cerr<<"First entry have no elements"<<std::endl;
    return false; //?
  }
  for(int i=0;i<tbl.size();++i)
  {
    if(tbl[i].size()!=sz)
    {
      std::cerr<<"Line :"<<i<<" have inconsistent number of entries:"
          <<tbl[i].size()<<" expected:"<<sz<<std::endl;
      for(int j=0;j<tbl[i].size();j++)
        std::cerr<<tbl[i][j].c_str()<<" - ";
      std::cerr<<std::endl;
      return false;
    }
  }
  return true;
}

class CommandProgressUpdate : public itk::Command
{
  public:
    typedef CommandProgressUpdate    Self;
    typedef itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer; 
    itkNewMacro( Self );
    std::vector<bool>               _progress;
  protected:
    CommandProgressUpdate():
      _progress(11,false)
    {};
    
  public:
    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      const itk::ProcessObject * filter =
          dynamic_cast< const itk::ProcessObject * >( object );
      if( ! itk::ProgressEvent().CheckEvent( &event ) )
      {
        return;
      }
      int prg=filter->GetProgress()*(_progress.size()-1);
      if(prg>=_progress.size()) prg=_progress.size()-1;
      
      if(!_progress[prg])
      {
        std::cout << prg*100.0/(_progress.size()-1) << "% " <<std::flush;
        _progress[prg]=true;
      } else {
        std::cout << "."<<std::flush;
      }
    }
    
    void ResetProgress(void)
    {
      _progress=std::vector<bool>(11,false);
      std::cout<<std::endl;
    }
};

void show_usage ( const char *name )
{
  std::cerr
      << "Usage: "<<name<<"  <input> <output normalized>" << std::endl
      << "\t--train <training_data> use this samples list" << std::endl
      << "\t--train2 <training_data2> use this samples list" << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--mask <mask>" <<std::endl
      << "\t--patch <n>" <<std::endl
      << "\t--search <n>" <<std::endl
      << "\t--beta <f> scaling coefficient for minimal distance, default 0.125 to be compatible with Pierric's code"<<std::endl
      << "\t--iter <n> maximum number of iterations when preselection is used"<<std::endl
      << "\t--extract <n> extract label N from the label image first"<<std::endl
      << "\t--top <n> select only top n samples from each training set"<<std::endl
      << "\t--float store results as floats"<<std::endl
      << "\t--short store results as short" <<std::endl
      << "\t--byte store results as  byte"  <<std::endl
      << "\t--box use Box structuring element (default)" <<std::endl
      << "\t--ball use ball structuring element "<<std::endl
      << "\t--median <n> median filter radius (default 4)" << std::endl;
}

//! a helper function for image saving
template<class T> void save_image(typename T::Pointer image,const char *fname,const char* minc_datatype,itk::Object* metadata,const std::string& history)
{
  if(metadata) itk::copy_metadata(image,metadata);
  if(!history.empty()) itk::append_minc_history(image,history);
  itk::set_minc_storage_type(image,minc_datatype);
  typename itk::ImageFileWriter< T >::Pointer writer = itk::ImageFileWriter<T>::New();
  writer->SetFileName(fname);
  writer->SetInput( image );
  writer->Update();
}

//! a helper function for image loading
template <class T> typename T::Pointer load_image(const char *file)
{
  typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();
  reader->SetFileName(file);
  reader->Update();
  return reader->GetOutput();
}

struct feature_image_distance
{
  FeatureImageType::Pointer img;
  double similarity;
  
  feature_image_distance()
  {}
  
  feature_image_distance(FeatureImageType::Pointer _img,double _dist=0):
    img(_img),similarity(_dist)
  {}
};

//! helper function that calculates SSD distances between all images in the train_images and in_img withing ROI
void ssd_distance(std::vector<feature_image_distance>& train_images,
                  FeatureImageType::Pointer in_img,
                  itk::mask3d::Pointer roi_in,
                  bool verbose)
{
  if(verbose)
    std::cout<<"Calculating SSD.."<<std::endl;
  
  for(size_t i=0;i<train_images.size();i++)
  {
    //let's do it old fashioned way, for now....
    itk::ImageRegionIterator<FeatureImageType> img_it(in_img, in_img->GetLargestPossibleRegion());
    itk::ImageRegionIterator<FeatureImageType> sample_it(train_images[i].img, in_img->GetLargestPossibleRegion());
    double dist=0.0;
    int cnt=0;
    
    img_it.GoToBegin();
    sample_it.GoToBegin();
    if(roi_in.IsNotNull())
    {
      itk::ImageRegionIterator<itk::mask3d> mask_it(roi_in, roi_in->GetLargestPossibleRegion());
      mask_it.GoToBegin();
      
      while(!img_it.IsAtEnd() && !sample_it.IsAtEnd() && !mask_it.IsAtEnd())
      {
        if(mask_it.Get()>0)
        {
          double d=(img_it.Get()-sample_it.Get());
          dist+=d*d;
          cnt++;
        }
        ++img_it;
        ++sample_it;
        ++mask_it;
      }
    } else {
      while(!img_it.IsAtEnd() && !sample_it.IsAtEnd())
      {
        double d=(img_it.Get()-sample_it.Get());
        dist+=d*d;
        cnt++;
        ++img_it;
        ++sample_it;
      }
    }
    if(cnt>0) dist/=cnt;
    //train_images[i].similarity=dist;
    train_images[i].similarity=dist;
  }
}

bool ssd_compare(const feature_image_distance &i,const feature_image_distance &j)
{
  return i.similarity < j.similarity;
}

int main ( int argc,char **argv )
{
  int clobber=0;
  int store_float=0;
  int store_short=0;
  int store_byte=0;
  int verbose=0;
  std::string train_f,train_f2;
  int selected;
  std::string mask_f;
  double sigma=0.0;
  int patch_radius=1;
  int search_radius=1;
  int select_top=0;
  int median_radius=4;
  
  int max_iterations=10;
  int ball_structuring_element=0;
  
  std::string history=itk::minc_time_stamp(argc, argv); 

  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose,1},
    {"quiet", no_argument, &verbose,  0},
    {"clobber", no_argument, &clobber,1},
    {"train", required_argument,      0, 't'},
    {"train2", required_argument,     0, 'T'},
    {"threshold", required_argument,  0, 'h'},
    {"mask", required_argument,       0, 'm'},
    {"patch", required_argument,      0, 'p'},
    {"search", required_argument,     0, 's'},
    {"beta", required_argument,       0, 'b'},
    {"iter", required_argument,       0, 'I'},
    {"median", required_argument,       0, 'M'},
    
    {"float",   no_argument, &store_float, 1},
    {"short",   no_argument, &store_short, 1},
    {"byte",    no_argument, &store_byte,  1},
    {"top", required_argument,        0, 'o'},
    
    {"box", no_argument, &ball_structuring_element, 0},
    {"ball", no_argument, &ball_structuring_element, 1},
    {0, 0, 0, 0}
  };

  int c;
  for ( ;; )
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long ( argc, argv, "t:r:m:T:S:C:G:I:M:", long_options, &option_index );

    /* Detect the end of the options. */
    if ( c == -1 )
      break;

    switch ( c )
    {
      case 0:
        break;
      case 't':
        train_f=optarg;
        break;
      case 'T':
        train_f2=optarg;
        break;
      case 'm':
        mask_f=optarg;
        break;
      case 'M':
        median_radius=atoi(optarg);
        break;
      case 'p':
        patch_radius=atoi(optarg);
        break;
      case 's':
        search_radius=atoi(optarg);
        break;
      case 'I':
        max_iterations=atoi(optarg);
        break;
      case 'o':
        select_top=atoi(optarg);
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage ( argv[0] );
        return 1;
    }
  }

  if ( ( argc - optind ) < 2 )
  {
    show_usage ( argv[0] );
    return 1;
  }

  std::string input_f=argv[optind];
  std::string output_f=argv[optind+1];

  strings inputs_f ( 1 );

  inputs_f[0]=input_f;

  if (!clobber  && !access (output_f.c_str(), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }

  try
  {
    FeatureImageType::Pointer in_img(FeatureImageType::New());
    
    in_img=load_image<image3d>(input_f.c_str());

    mask3d::Pointer roi_in(mask3d::New());
    
    if(!mask_f.empty())
      roi_in=load_image<itk::mask3d>(mask_f.c_str());
    
    string_table tbl;
    string_table tbl2;
    
    read_table_n ( train_f.c_str(),tbl,2 );
    if ( !table_verify ( tbl ) )
      return 1;
    
    //TODO: figure out if we need to adjust weights for the addditional samples 
    if(!train_f2.empty())
    {
      if(verbose)
        std::cout<<"Loading secondary training set!"<<std::endl;
      
      read_table_n ( train_f2.c_str(),tbl2,2 );
      
      if ( !table_verify ( tbl2 ) )
        return 1;
    }
    
    if(verbose)
      std::cout<<"Loading volumes..."<<std::endl;

    std::vector<feature_image_distance> train_images,train_images2;
    std::vector<FeatureImageType::Pointer> trg_images;

    for(size_t i=0;i<tbl.size();i++)
    {

      if(verbose) 
        std::cout<<tbl[i][0].c_str()<<"\t"<<std::flush;
      FeatureImageType::Pointer feature=load_image<FeatureImageType>(tbl[i][0].c_str());

      if(verbose) 
        std::cout<<tbl[i][1].c_str()<<"\t"<<std::flush;

      if(verbose) 
        std::cout<<std::endl;

      train_images.push_back(feature_image_distance(feature));
    }

    for(size_t i=0;i<tbl2.size();i++)
    {

      if(verbose) 
        std::cout<<tbl2[i][0].c_str()<<"\t"<<std::flush;
      FeatureImageType::Pointer feature=load_image<FeatureImageType>(tbl2[i][0].c_str());

      if(verbose) 
        std::cout<<tbl2[i][1].c_str()<<"\t"<<std::flush;

      if(verbose) 
        std::cout<<std::endl;

      train_images2.push_back(feature_image_distance(feature));
    }
    
    
    if(select_top)
    {
      // calculate SSD similarities
      ssd_distance(train_images ,in_img,roi_in,verbose);
      ssd_distance(train_images2,in_img,roi_in,verbose);
      
      std::sort(train_images.begin(), train_images.end(), ssd_compare);
      std::sort(train_images2.begin(),train_images2.end(),ssd_compare);
      
      if(train_images.size()>select_top)
        train_images.erase(train_images.begin()+select_top,train_images.end()); 

      if(train_images2.size()>select_top)
        train_images2.erase(train_images2.begin()+select_top,train_images2.end()); 
    }

    //unite training sets
    for(size_t i=0;i<train_images2.size();i++)
      trg_images.push_back(train_images2[i].img);

    //hope to reclaim some memory (?)
    train_images2.clear();

    //put images into preselection list
    for(size_t i=0;i<train_images.size();i++)
      trg_images.push_back(train_images[i].img);

    if(verbose)
      std::cout<<"done"<<std::endl;

    //2. calculate local similarities
    if(verbose)
      std::cout<<"Patch based analysis, "<<std::endl
              <<"\tpatch size="<<patch_radius<<" "<<std::endl
              <<"\tsearch radius="<<search_radius<<" "<<std::endl
              <<"\tsamples="<<trg_images.size()<<" "<<std::endl
              <<"\tUsing "<<(ball_structuring_element?"Ball":"Box")<<" structural element"<<std::endl
              <<"\tMedian Filter radius "<<median_radius<<std::endl;

    StructuringElementType::RadiusType search_radius_={{search_radius,search_radius,search_radius}};
    StructuringElementType::RadiusType patch_radius_={{patch_radius,patch_radius,patch_radius}};

    StructuringElementType  searchKernel=ball_structuring_element?StructuringElementType::Box(search_radius_):StructuringElementType::Ball(search_radius_);
    StructuringElementType  patchKernel=ball_structuring_element?StructuringElementType::Box(patch_radius_):StructuringElementType::Ball(patch_radius_);

    NormalizationType::Pointer normalization=NormalizationType::New();
    MedianFilterType::Pointer  median_filter=MedianFilterType::New();
    
    MedianFilterType::RadiusType median_radius_sz;
    median_radius_sz.Fill(median_radius);
    
    median_filter->SetRadius(median_radius_sz);

    NormalizationFilterNoPsType::Pointer filter(NormalizationFilterNoPsType::New());
    
    FeatureImageType::Pointer corr_img(FeatureImageType::New());
    allocate_same<FeatureImageType,FeatureImageType>(corr_img,in_img);
    corr_img->FillBuffer(0.0);
    
    NegExpFilterType::Pointer neg_exp_filter=NegExpFilterType::New();
    MultiplyFilterType::Pointer multiply_filter=MultiplyFilterType::New();
    AddFilterType::Pointer add_filter=AddFilterType::New();
    typedef itk::StatisticsImageFilter<FeatureImageType> StatisticsImageFilterType;
    StatisticsImageFilterType::Pointer stat_filter=StatisticsImageFilterType::New();
    
    
    filter->SetNormalizationLibrary(trg_images);
    filter->SetSearchKernel(searchKernel);
    filter->SetPatchKernel(patchKernel);
    filter->SetProcess(normalization);
    
    if(!mask_f.empty())
      filter->SetRoiImage(roi_in);

    CommandProgressUpdate::Pointer observer = CommandProgressUpdate::New();
    filter->AddObserver(itk::ProgressEvent(),observer);

    stat_filter->SetInput(corr_img);
    stat_filter->Update();

    for(int it=0;it<max_iterations;it++)
    {
      if(verbose)
        std::cout<<"Iteration "<<it<<" of "<<max_iterations<<std::endl;
      
      corr_img->Modified();
      in_img->Modified();
      
      neg_exp_filter->SetInput(corr_img);
      multiply_filter->SetInput1(neg_exp_filter->GetOutput());
      multiply_filter->SetInput2(in_img);

      filter->SetInput(multiply_filter->GetOutput());
      median_filter->SetInput(filter->GetOutput());
      
      add_filter->SetInput1(corr_img);
      add_filter->SetInput2(median_filter->GetOutput());

      add_filter->Update();
      
      corr_img=add_filter->GetOutput();
      
      stat_filter->SetInput(median_filter->GetOutput());
      stat_filter->Update();
      
      if(verbose)
      {
        observer->ResetProgress();
        std::cout<<" Update mean:"<< stat_filter->GetMean() << std::endl;
        std::cout<<" Update std:"<<  stat_filter->GetSigma() << std::endl;
      }
    }

    if(verbose)
      std::cout<<"done"<<std::endl;

    save_image<FeatureImageType>(corr_img,output_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);

  } catch( itk::ExceptionObject & err ) 
  {
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return 2;
  } 
  return 0;
}

// kate: space-indent on; hl c++;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 2 
