/* ----------------------------- MNI Header -----------------------------------
@NAME       :  itk_patch_morphology
@DESCRIPTION:  patch morphology
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
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>

#include "itkSegmentationNonLocalFilter.h"
#include "itkSegmentationPreselectionFilter.h"

#include <iostream>
#include <fstream>
#include <set>
#include <getopt.h>
#include <libgen.h>
#include <unistd.h>

typedef itk::Image<float,3> FeatureImageType;
typedef itk::Image<unsigned char,3> LabelImageType;

typedef itk::ImagePair<FeatureImageType,FeatureImageType> ImagePairType;
typedef itk::FlatStructuringElement< 3  >            StructuringElementType;

typedef itk::MinWeightedAverageLabeling<float,float,itk::InvExpWeight<double>,3 >          WeightedAverageLabelingType;
typedef itk::MinWeightedDiscreteLabeling<float,unsigned char,itk::InvExpWeight<double>,3 > DiscreteLabelingType;
//typedef itk::MedianWeightedAverageLabeling<float,float,itk::InvExpWeight<double>,3 > MedianWeightedAverageLabelingType;
//typedef itk::MinimalDistanceNonLocalFilter<float,float,itk::InvExpWeight<double>,3 > MinWeightedAverageLabelingType;


typedef itk::SegmentationMeanAndSdPreselectionFilter<FeatureImageType,StructuringElementType>  MeanAndSdPreselectionType;

typedef itk::SegmentationNonLocalFilter<
                            FeatureImageType,
                            FeatureImageType,
                            FeatureImageType,
                            StructuringElementType,
                            StructuringElementType,
                            WeightedAverageLabelingType,
                            itk::L2PatchDistance<FeatureImageType,StructuringElementType>,
                            MeanAndSdPreselectionType>  
                            SegmentationFilterType;

typedef itk::SegmentationNonLocalFilter<
                            FeatureImageType,
                            FeatureImageType,
                            itk::mask3d,
                            StructuringElementType,
                            StructuringElementType,
                            DiscreteLabelingType,
                            itk::L2PatchDistance<FeatureImageType,StructuringElementType>,
                            MeanAndSdPreselectionType>
                            DiscreteSegmentationFilterType;

typedef itk::SegmentationNonLocalFilter<
                            FeatureImageType,
                            FeatureImageType,
                            FeatureImageType,
                            StructuringElementType,
                            StructuringElementType,
                            WeightedAverageLabelingType,
                            itk::L2PatchDistance<FeatureImageType,StructuringElementType> >  
                            SegmentationFilterNoPsType;

typedef itk::SegmentationNonLocalFilter<
                            FeatureImageType,
                            FeatureImageType,
                            itk::mask3d,
                            StructuringElementType,
                            StructuringElementType,
                            DiscreteLabelingType,
                            itk::L2PatchDistance<FeatureImageType,StructuringElementType> >
                            DiscreteSegmentationFilterNoPsType;

typedef itk::BinaryThresholdImageFilter<FeatureImageType,FeatureImageType> LabelThresholdFilter;
typedef itk::BinaryThresholdImageFilter<itk::mask3d,FeatureImageType> LabelThresholdFilterD;

typedef SegmentationFilterType::Superclass BaseFilterType;
typedef DiscreteSegmentationFilterType::Superclass DiscreteBaseFilterType;

using namespace std;

//! vector of strings 
typedef std::vector<std::string> strings;

//! string table (vector of vectors of strings)
typedef std::vector<strings> string_table;

FeatureImageType::Pointer extract_label_in_image(FeatureImageType::Pointer img,int label_to_extract)
{
  LabelThresholdFilter::Pointer threshold=LabelThresholdFilter::New();
  threshold->SetInPlace(true);
  threshold->SetLowerThreshold(label_to_extract-0.5);
  threshold->SetUpperThreshold(label_to_extract+0.5);
  threshold->SetInsideValue(1.0);
  threshold->SetOutsideValue(0.0);
  threshold->SetInput(img);
  threshold->Update();
  return threshold->GetOutput();
}

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
    typedef CommandProgressUpdate   Self;
    typedef itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer; 
    itkNewMacro( Self );
    std::vector<bool> 	  						  _progress;
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

bool ssd_compare(const ImagePairType &i,const ImagePairType &j)
{
  return i.similarity < j.similarity;
}

void show_usage ( const char *name )
{
  std::cerr
      << "Patch-based segmentation tool" << std::endl
      << "Reference: Patch-based segmentation using expert priors: Application to hippocampus and ventricle segmentation" << std::endl
      << "\thttp://dx.doi.org/10.1016/j.neuroimage.2010.09.018" << std::endl << std::endl
      << "Usage: "<<name<<"  <input> [output labels]" << std::endl
      << "\t--train <training_data> use this samples list" << std::endl
      << "\t--train2 <training_data2> use this samples list" << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--mask <mask>" <<std::endl
      << "\t--patch <n>" <<std::endl
      << "\t--search <n>" <<std::endl
      //<< "\t--dist <output> output minimal patch distance" <<std::endl
      << "\t--cls  <output> output closest patch id" <<std::endl
      << "\t--beta <f> scaling coefficient for minimal distance, default 0.125 to be compatible with Pierric's code"<<std::endl
      << "\t--scaling <f> use this map for scaling"<<std::endl
      << "\t--threshold <d> preselection threshold 0 - select everything, 1.0 - select nothing"<<std::endl
      << "\t--discrete <n> use for discrete labeling with n number of classes, including  0"<<std::endl
      << "\t--confidence <output> output confidence map"<<std::endl
      << "\t--adist <output> output mean distance map"<<std::endl
      << "\t--grading <output> output image grading map"<<std::endl
      << "\t--iter <n> maximum number of iterations when preselection is used"<<std::endl
      << "\t--extract <n> extract label N from the label image first"<<std::endl
      << "\t--top <n> select only top n samples from each training set"<<std::endl
      << "\t--float store results as floats"<<std::endl
      << "\t--short store results as short" <<std::endl
      << "\t--byte store results as  byte"  <<std::endl
      << "\t--box use Box structuring element (default)" <<std::endl
      << "\t--ball use ball structuring element "<<std::endl
      << "\t--prelabel <labels> prelabel data using this labels, segment only parts which are not labeled yet"<<std::endl
      << "\t--prob <prefix> store per-label probabilities"<<std::endl
      << "\t--groups <n> - modifies behaviour to allow one more column in the training sample (group id) and gurantees balanced sample after pre-selection"<<std::endl;
}

//! a helper function for image saving
template<class T> void save_image(typename T::Pointer image,const char *fname,const char* minc_datatype,itk::Object* metadata,const std::string& history)
{
  if(metadata) copy_metadata(image,metadata);
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

//! helper function that calculates SSD distances between all images in the train_images and in_img withing ROI
void ssd_distance(std::vector<ImagePairType>& train_images ,FeatureImageType::Pointer in_img,itk::mask3d::Pointer roi_in,bool verbose)
{
  if(verbose)
    std::cout<<"Calculating SSD.."<<std::endl;
  
  for(size_t i=0;i<train_images.size();i++)
  {
    //let's do it old fashioned way, for now....
    itk::ImageRegionIterator<FeatureImageType> img_it(in_img, in_img->GetLargestPossibleRegion());
    itk::ImageRegionIterator<FeatureImageType> sample_it(train_images[i].feature, in_img->GetLargestPossibleRegion());
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
    train_images[i].similarity=dist;
  }
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
  double beta=0.125;
  int extract_label=0;
  int select_top=0;
  std::string output_dist_f,output_cls_f,output_conf_f,output_adist_f,output_grading_f;
  std::string prob_prefix_f;
  std::string preload_f;
  std::string history=itk::minc_time_stamp(argc, argv); 
  double preselect_threshold=0.0;
  int discrete=0;
  int max_iterations=50;
  int ball_structuring_element=0;
  int groups=0;

  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"train", required_argument, 0, 't'},
    {"train2", required_argument, 0, 'T'},
    {"threshold", required_argument, 0, 'h'},
    {"mask", required_argument, 0, 'm'},
    {"patch", required_argument, 0, 'p'},
    {"search", required_argument, 0, 's'},
    {"beta", required_argument, 0, 'b'},
    //{"dist", required_argument, 0, 'd'},
    {"discrete", required_argument, 0, 'D'},
    {"extract",required_argument, 0, 'E'},
    {"iter", required_argument, 0, 'I'},
    
    {"cls", required_argument, 0, 'c'},
    {"roi",     required_argument, 0,    'r'},
    {"float",   no_argument, &store_float, 1},
    {"short",   no_argument, &store_short, 1},
    {"byte",    no_argument, &store_byte,  1},
    {"confidence", required_argument, 0, 'C'},
    {"adist", required_argument, 0, 'S'},
    {"grading", required_argument, 0, 'G'},
    {"top", required_argument, 0, 'o'},
    {"prelabel", required_argument, 0, 'P'},
    {"prob", required_argument, 0, 'R'},
    {"groups", required_argument, 0, 'g'},
    
    {"box", no_argument, &ball_structuring_element, 0},
    {"ball", no_argument, &ball_structuring_element, 1},
    {0, 0, 0, 0}
  };

  int c;
  for ( ;; )
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long ( argc, argv, "t:r:m:T:D:S:C:G:I:E:P:", long_options, &option_index );

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
      case 'p':
        patch_radius=atoi(optarg);
        break;
      case 's':
        search_radius=atoi(optarg);
        break;
      case 'D':
        discrete=atoi(optarg);
        break;
      case 'I':
        max_iterations=atoi(optarg);
        break;
      case 'b':
        beta=atof(optarg);
        break;
//       case 'd':
//         output_dist_f=optarg;
//         break;
      case 'c':
        output_cls_f=optarg;
        break;
      case 'g':
        groups=atoi(optarg);
        break;
      case 'h':
        preselect_threshold=atof(optarg);
        break;
      case 'C':
        output_conf_f=optarg;
        break;
      case 'S':
        output_adist_f=optarg;
        break;
      case 'G':
        output_grading_f=optarg;
        break;
      case 'E':
        extract_label=atoi(optarg);
        break;
      case 'o':
        select_top=atoi(optarg);
        break;
      case 'P':
        preload_f=optarg;
        break;
      case 'R':
        prob_prefix_f=optarg;
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage ( argv[0] );
        return 1;
    }
  }

  if ( ( argc - optind ) < 1 )
  {
    show_usage ( argv[0] );
    return 1;
  }

  std::string input_f=argv[optind];
  std::string output_f;

  if(( argc - optind ) > 1)
    output_f=argv[optind+1];


  strings inputs_f ( 1 );

  inputs_f[0]=input_f;

  if (!clobber && !output_f.empty() && !access (output_f.c_str(), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }

  if (!clobber && !output_dist_f.empty() && !access (output_dist_f.c_str(), F_OK))
  {
    std::cerr << output_dist_f.c_str () << " Exists!" << std::endl;
    return 1;
  }

  if (!clobber && !output_cls_f.empty() && !access (output_cls_f.c_str(), F_OK))
  {
    std::cerr << output_cls_f.c_str () << " Exists!" << std::endl;
    return 1;
  }
  
  if (!clobber && !output_conf_f.empty() && !access (output_conf_f.c_str(), F_OK))
  {
    std::cerr << output_conf_f.c_str () << " Exists!" << std::endl;
    return 1;
  }
  
  if (!clobber && !output_adist_f.empty() && !access (output_adist_f.c_str(), F_OK))
  {
    std::cerr << output_adist_f.c_str () << " Exists!" << std::endl;
    return 1;
  }
  
  if (!clobber && !output_grading_f.empty() && !access (output_grading_f.c_str(), F_OK) && !groups)
  {
    std::cerr << output_grading_f.c_str () << " Exists!" << std::endl;
    return 1;
  }
    
  try
  {
    FeatureImageType::Pointer in_img(FeatureImageType::New());
    FeatureImageType::Pointer out_conf,out_adist,out_grading;
    FeatureImageType::Pointer preload,preload_confidence;
    itk::mask3d::Pointer   preload_discrete;
    
    
    in_img=load_image<itk::image3d>(input_f.c_str());

    itk::mask3d::Pointer roi_in(itk::mask3d::New());
    
    if(!mask_f.empty())
      roi_in=load_image<itk::mask3d>(mask_f.c_str());
    
    if(!preload_f.empty())
    {
      
    if(!discrete) //select anything above 0.5 to have 1.0 confidence
    {
        preload=load_image<itk::image3d>(preload_f.c_str());
        LabelThresholdFilter::Pointer threshold=LabelThresholdFilter::New();
        threshold->SetInPlace(false);
        threshold->SetLowerThreshold(0.5);
        threshold->SetUpperThreshold(1.5);
        threshold->SetInsideValue(1.0);
        threshold->SetOutsideValue(0.0);
        threshold->SetInput(preload);
        threshold->Update();
        preload_confidence=threshold->GetOutput();
        preload_confidence->DisconnectPipeline();
    } else {
        preload_discrete=load_image<itk::mask3d>(preload_f.c_str());
        LabelThresholdFilterD::Pointer threshold=LabelThresholdFilterD::New();
        threshold->SetInPlace(false);
        threshold->SetLowerThreshold(1);
        threshold->SetUpperThreshold(discrete);
        threshold->SetInsideValue(1.0);
        threshold->SetOutsideValue(0.0);
        threshold->SetInput(preload_discrete);
        threshold->Update();
        preload_confidence=threshold->GetOutput();
        preload_confidence->DisconnectPipeline();
    }
    }
    
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
    
    std::vector<ImagePairType> train_images,train_images2;
    std::vector<FeatureImageType::Pointer> trg_images;
    
    if(!groups)
    {
      for(size_t i=0;i<tbl.size();i++)
      {
          
        if(verbose) 
          std::cout<<tbl[i][0].c_str()<<"\t"<<std::flush;
        FeatureImageType::Pointer feature=load_image<FeatureImageType>(tbl[i][0].c_str());
        
        if(verbose) 
          std::cout<<tbl[i][1].c_str()<<"\t"<<std::flush;
        
        FeatureImageType::Pointer label=load_image<FeatureImageType>(tbl[i][1].c_str());
        
        double grading=1;
        
        if(tbl[i].size()>2)
        {
          grading=atof(tbl[i][2].c_str());
          if(verbose)
            std::cout<<"\t"<<tbl[i][2].c_str();
        }
        if(verbose) 
          std::cout<<std::endl;
          
        if(extract_label>0 ) //extract label here
        {
          FeatureImageType::Pointer label2=extract_label_in_image(label,extract_label);
          train_images.push_back(ImagePairType(feature,label2,grading));
        } else 
          train_images.push_back(ImagePairType(feature,label,grading));
      }
      
      for(size_t i=0;i<tbl2.size();i++)
      {
        if(verbose) 
          std::cout<<tbl2[i][0].c_str()<<"\t"<<std::flush;
        FeatureImageType::Pointer feature=load_image<FeatureImageType>(tbl2[i][0].c_str());
        
        if(verbose) 
          std::cout<<tbl2[i][1].c_str()<<"\t"<<std::flush;
        
        FeatureImageType::Pointer label=load_image<FeatureImageType>(tbl2[i][1].c_str());
        
        double grading=1.0;
        
        if(tbl2[i].size()>2)
        {
          grading=atof(tbl2[i][2].c_str());
          if(verbose)
            std::cout<<"\t"<<tbl2[i][2].c_str();
        }
        if(verbose) 
          std::cout<<std::endl;
          
        if(extract_label>0 ) //extract label here
        {
          FeatureImageType::Pointer label2=extract_label_in_image(label,extract_label);
          train_images2.push_back(ImagePairType(feature,label2,grading));
        } else 
          train_images2.push_back(ImagePairType(feature,label,grading));
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
        train_images.push_back(train_images2[i]);
      
      //hope to reclaim some memory (?)
      train_images2.clear();
    }  else {
      std::set<int> group_set;
      for(size_t i=0;i<tbl.size();i++)
      {
        if(verbose) 
          std::cout<<tbl[i][0].c_str()<<"\t"<<std::flush;
          
        FeatureImageType::Pointer feature=load_image<FeatureImageType>(tbl[i][0].c_str());
        
        if(verbose) 
          std::cout<<tbl[i][1].c_str()<<"\t"<<std::flush;
        
        FeatureImageType::Pointer label=load_image<FeatureImageType>(tbl[i][1].c_str());
        
        double grading=1;
        int    group=-1;
        
        if(tbl[i].size()>2)
        {
          grading=atof(tbl[i][2].c_str());
          if(verbose)
            std::cout<<"\t"<<grading;
        }
        if(verbose) 
          std::cout<<std::endl;
          
        if(tbl[i].size()>3)
        {
          group=atoi(tbl[i][3].c_str());
          if(verbose)
            std::cout<<"\t"<<group;
        }
        
        if(verbose) 
          std::cout<<std::endl;
        
        group_set.insert(group);
        
        if(extract_label>0 ) //extract label here
        {
          FeatureImageType::Pointer label2=extract_label_in_image(label,extract_label);
          train_images.push_back(ImagePairType(feature,label2,grading,0,group));
        } else 
          train_images.push_back(ImagePairType(feature,label,grading,0,group));
      }
      if(select_top)
      {
        // calculate SSD similarities
        ssd_distance(train_images,in_img,roi_in,verbose);
        std::vector<ImagePairType> _train_images;
        for(std::set<int>::iterator i=group_set.begin();i!=group_set.end();++i)
        {
          std::vector<ImagePairType> grp_images;
          for(std::vector<ImagePairType>::iterator j=train_images.begin();j!=train_images.end();++j)
          {
            if((*j).group==(*i)) 
              grp_images.push_back((*j));
          }
          std::sort(grp_images.begin(), grp_images.end(), ssd_compare);
          for(int k=0;k<select_top && k< grp_images.size(); k++)
            _train_images.push_back(grp_images[k]);
        }
        train_images=_train_images;
      }
    }

    //put images into preselection list
    for(size_t i=0;i<train_images.size();i++)
      trg_images.push_back(train_images[i].feature);

    if(verbose)
      std::cout<<"done"<<std::endl;
    

    //2. calculate local similarities
    if(verbose)
      std::cout<<"Patch based analysis, "<<std::endl
              <<"\tpatch size="<<patch_radius<<" "<<std::endl
              <<"\tsearch radius="<<search_radius<<" "<<std::endl
              <<"\tsamples="<<train_images.size()<<" "<<std::endl
              <<"\tpreselection threshold="<<preselect_threshold<<std::endl
              <<"\tSmoothing parameter Beta="<<beta<<std::endl
              <<"\tUsing "<<(ball_structuring_element?"Ball":"Box")<<" structural element"<<std::endl;


    MeanAndSdPreselectionType::Pointer   preselection=MeanAndSdPreselectionType::New();
    StructuringElementType::RadiusType search_radius_;search_radius_.assign(search_radius);//={{search_radius,search_radius,search_radius}};
    StructuringElementType::RadiusType patch_radius_;patch_radius_.assign(patch_radius);//={{patch_radius,patch_radius,patch_radius}};

    StructuringElementType  searchKernel=ball_structuring_element?StructuringElementType::Box(search_radius_):StructuringElementType::Ball(search_radius_);
    StructuringElementType  patchKernel=ball_structuring_element?StructuringElementType::Box(patch_radius_):StructuringElementType::Ball(patch_radius_);

    WeightedAverageLabelingType::Pointer weightedAverage=WeightedAverageLabelingType::New();
    DiscreteLabelingType::Pointer discreteWeightedAverage=DiscreteLabelingType::New();

    weightedAverage->GetWeight().SetBeta(beta); // finally, compatible with Pierrick's code
    discreteWeightedAverage->GetWeight().SetBeta(beta); // finally, compatible with Pierrick's code

    if(!discrete) //we are dealing with continious labeling (i.e fuzzy segmentation case)
    {
      BaseFilterType::Pointer flt;
      CommandProgressUpdate::Pointer observer = CommandProgressUpdate::New();
      if(preselect_threshold>0) // we have preselection threshold
      {
        SegmentationFilterType::Pointer filter(SegmentationFilterType::New());
        
        if(verbose)
          std::cout<<"Calculating preselection information..."<<std::flush;
        
        preselection->SetThreshold(preselect_threshold);
        preselection->SetPatchKernel(patchKernel);
        preselection->SetImages(in_img,trg_images);
        
        if(verbose)
          std::cout<<"done"<<std::endl;
        
        filter->SetPreselectionFilter(preselection);
        filter->SetSegmentationLibrary(train_images);
        
        filter->SetPatchKernel(patchKernel);
        filter->SetProcess(weightedAverage);
        
        if(!mask_f.empty())
          filter->SetRoiImage(roi_in);
        
        if(!preload_f.empty())
        {
          filter->SetPreloadedResults(preload);
          filter->SetConfidence(preload_confidence);
        }
        
        
        filter->SetInput(in_img);
        filter->AddObserver(itk::ProgressEvent(),observer);
        
        size_t iter=0;
        do
        {
          filter->SetSearchKernel(searchKernel);
          filter->Update();
          
          iter++;
          observer->ResetProgress();
          if(verbose)
            std::cout<<"Iteration:"<<iter<<" number of non-confident voxels:"<<filter->GetNonConfidentVoxels()<<" at threshold:"<<preselect_threshold<<" Search radius:"<<search_radius<<std::endl;
          
          preselect_threshold*=0.95;
          preselection->SetThreshold(preselect_threshold);
          filter->SetPreselectionFilter(preselection);
          
          search_radius++;
          
          StructuringElementType::RadiusType search_radius_;search_radius_.assign(search_radius);//={{search_radius,search_radius,search_radius}};
          searchKernel=ball_structuring_element?StructuringElementType::Box(search_radius_):StructuringElementType::Ball(search_radius_);
        } while(iter<max_iterations && filter->GetNonConfidentVoxels()>0);

        out_conf=filter->GetConfidence();
        out_adist=filter->GetSearchDistance();
        out_grading=filter->GetGradingMap();
        
        flt=filter;
      } else {
        SegmentationFilterNoPsType::Pointer filter(SegmentationFilterNoPsType::New());
        
        filter->SetSegmentationLibrary(train_images);
        filter->SetSearchKernel(searchKernel);
        filter->SetPatchKernel(patchKernel);
        filter->SetProcess(weightedAverage);
        if(!mask_f.empty())
          filter->SetRoiImage(roi_in);
        
        if(!preload_f.empty())
        {
          filter->SetPreloadedResults(preload);
          filter->SetConfidence(preload_confidence);
        }
        
        filter->SetInput(in_img);
        filter->AddObserver(itk::ProgressEvent(),observer);
        filter->Update();
        
        out_conf=filter->GetConfidence();
        out_adist=filter->GetSearchDistance();
        out_grading=filter->GetGradingMap();
        
        flt=filter;
      }
        
      if(verbose)
        std::cout<<"done"<<std::endl;

      if(!output_f.empty())
        save_image<FeatureImageType>(flt->GetOutput(),output_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);
      
      if(!output_conf_f.empty())
        save_image<FeatureImageType>(out_conf,output_conf_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);
      
      if(!output_adist_f.empty())
        save_image<FeatureImageType>(out_adist,output_adist_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);
      
      if(!output_grading_f.empty())
        save_image<FeatureImageType>(out_grading,output_grading_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);

    } else { //USING discrete labeling mode
      discreteWeightedAverage->SetLabelCount(discrete);
      
      DiscreteBaseFilterType::Pointer flt;
      CommandProgressUpdate::Pointer observer = CommandProgressUpdate::New();
      if(preselect_threshold>0)
      {
        DiscreteSegmentationFilterType::Pointer filter(DiscreteSegmentationFilterType::New());
        
        if(verbose)
          std::cout<<"Calculating preselection information..."<<std::flush;
        
        preselection->SetThreshold(preselect_threshold);
        preselection->SetPatchKernel(patchKernel);
        preselection->SetImages(in_img,trg_images);
        
        if(verbose)
          std::cout<<"done"<<std::endl;
        
        filter->SetPreselectionFilter(preselection);
        filter->SetSegmentationLibrary(train_images);
        filter->SetPatchKernel(patchKernel);
        filter->SetProcess(discreteWeightedAverage);
        
        if(!mask_f.empty())
          filter->SetRoiImage(roi_in);
        
        if(!preload_f.empty())
        {
          filter->SetPreloadedResults(preload_discrete);
          filter->SetConfidence(preload_confidence);
        }
        
        filter->SetInput(in_img);
        filter->AddObserver(itk::ProgressEvent(),observer);
        
        size_t iter=0;
        do
        {
          filter->SetSearchKernel(searchKernel);
          filter->Update();
          
          iter++;
          observer->ResetProgress();
          
          if(verbose)
            std::cout<<"Iteration:"<<iter<<" number of non-confident voxels:"<<filter->GetNonConfidentVoxels()<<" at threshold:"<<preselect_threshold<<" Search radius:"<<search_radius<<std::endl;
          
          preselect_threshold*=0.95;
          preselection->SetThreshold(preselect_threshold);
          filter->SetPreselectionFilter(preselection);
          
          search_radius++;
          
          StructuringElementType::RadiusType search_radius_;search_radius_.assign(search_radius);//={search_radius,search_radius,search_radius};
          searchKernel=ball_structuring_element?StructuringElementType::Box(search_radius_):StructuringElementType::Ball(search_radius_);
          
        } while(iter<max_iterations && filter->GetNonConfidentVoxels()>0);
        
        out_conf=filter->GetConfidence();
        out_adist=filter->GetSearchDistance();
        out_grading=filter->GetGradingMap();
        
        if(!prob_prefix_f.empty())
        {
          for(int l=0;l<discrete;l++)
          {
            char tmp[1024];
            sprintf(tmp,"%s_%02d.mnc",prob_prefix_f.c_str(),l);
            save_image<FeatureImageType>(filter->GetProbabilityMap(l),tmp,store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);
          }
        }
        
        flt=filter;
        
      } else {
        DiscreteSegmentationFilterNoPsType::Pointer filter(DiscreteSegmentationFilterNoPsType::New());
        
        filter->SetSegmentationLibrary(train_images);
        filter->SetSearchKernel(searchKernel);
        filter->SetPatchKernel(patchKernel);
        filter->SetProcess(discreteWeightedAverage);
        
        if(!mask_f.empty())
          filter->SetRoiImage(roi_in);
        
        if(!preload_f.empty())
        {
          filter->SetPreloadedResults(preload_discrete);
          filter->SetConfidence(preload_confidence);
        }

        filter->SetInput(in_img);
        filter->AddObserver(itk::ProgressEvent(),observer);
        filter->Update();
        
        out_conf=filter->GetConfidence();
        out_adist=filter->GetSearchDistance();
        out_grading=filter->GetGradingMap();
        
        if(!prob_prefix_f.empty())
        {
          for(int l=0;l<discrete;l++)
          {
            char tmp[1024];
            sprintf(tmp,"%s_%02d.mnc",prob_prefix_f.c_str(),l);
            save_image<FeatureImageType>(filter->GetProbabilityMap(l),tmp,store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);            
          }
        }
        
        flt=filter;
      }
      
      if(verbose)
        std::cout<<"done"<<std::endl;
      if(!output_f.empty())
        save_image<itk::mask3d>(flt->GetOutput(),output_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);
      
      if(!output_conf_f.empty())
        save_image<FeatureImageType>(out_conf,output_conf_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);
      
      if(!output_adist_f.empty())
        save_image<FeatureImageType>(out_adist,output_adist_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);
      
      if(!output_grading_f.empty())
        save_image<FeatureImageType>(out_grading,output_grading_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);

    }
  } catch( itk::ExceptionObject & err ) 
  {
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return 2;
  } 
  return 0;
}

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2
