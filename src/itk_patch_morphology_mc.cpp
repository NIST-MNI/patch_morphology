/* ----------------------------- MNI Header -----------------------------------
@NAME       :  itk_patch_morphology_mc
@DESCRIPTION:  multi-channel patch morphology
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
#include <itkComposeImageFilter.h>

#include "itkSegmentationProcess.h"
#include "itkSegmentationNonLocalFilter.h"
#include "itkSegmentationPreselectionFilter.h"

#include <set>
#include <iostream>
#include <getopt.h>
#include <libgen.h>
#include <unistd.h>

typedef itk::Image<float,3>         SingleFeatureImageType;
typedef itk::VectorImage<float, 3 > MultipleFeatureImageType;
typedef itk::Image<unsigned char,3> LabelImageType;

typedef itk::ImagePair<MultipleFeatureImageType,SingleFeatureImageType> ImagePairType;
typedef itk::FlatStructuringElement< 3  >            StructuringElementType;

typedef itk::MinWeightedDiscreteLabeling<
              float,
              unsigned char,
              itk::InvExpWeight<double>,3 > 
              DiscreteLabelingType;

typedef itk::VectorL2PatchDistance< MultipleFeatureImageType,
                                    StructuringElementType> WeightedVectorPatchDistance;
                
typedef itk::SegmentationNonLocalFilter<
              MultipleFeatureImageType,
              SingleFeatureImageType,
              LabelImageType,
              StructuringElementType,
              StructuringElementType,
              DiscreteLabelingType,
              WeightedVectorPatchDistance >
                DiscreteSegmentationFilterNoPsType;

typedef itk::BinaryThresholdImageFilter<
              SingleFeatureImageType,
              SingleFeatureImageType > 
                LabelThresholdFilter;

typedef itk::BinaryThresholdImageFilter<
              LabelImageType,
              SingleFeatureImageType> 
                LabelThresholdFilterD;

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


#include <string>
#include <cstring>    // for strchr


/*****************************************************************
* This is the only part of the implementation that I don't like.
* It can probably be improved upon by the reader...
*/
namespace {
  inline bool
      isws (char c, char const * const wstr)
  {
    return (strchr(wstr,c) != NULL);
  }
}


/*****************************************************************
* Simplistic and quite Standard, but a bit slow.  This should be
* templatized on basic_string instead, or on a more generic StringT
* that just happens to support ::size_type, .substr(), and so on.
* I had hoped that "whitespace" would be a trait, but it isn't, so
* the user must supply it.  Enh, this lets them break up strings on
* different things easier than traits would anyhow.
*/
template <typename Container>
    void
    stringtok (Container &l, std::string const &s, char const * const ws = " \t\n")
{
  const std::string::size_type  S = s.size();
  std::string::size_type  i = 0;

  while (i < S) {
        // eat leading whitespace
    while ((i < S) && (isws(s[i],ws)))  ++i;
    if (i == S)  return;  // nothing left but WS

        // find end of word
    std::string::size_type  j = i+1;
    while ((j < S) && (!isws(s[j],ws)))  ++j;

        // add word
    l.push_back(s.substr(i,j-i));

        // set up for next loop
    i = j+1;
  }
}

template <typename Container>
    void
    stringtokd (Container &l, std::string const &s, char const * const ws = " \t\n")
{
  const std::string::size_type  S = s.size();
  std::string::size_type  i = 0;

  while (i < S) {
        // eat leading whitespace
    while ((i < S) && (isws(s[i],ws)))  ++i;
    if (i == S)  return;  // nothing left but WS

        // find end of word
    std::string::size_type  j = i+1;
    while ((j < S) && (!isws(s[j],ws)))  ++j;

        // add word
    l.push_back(atof(s.substr(i,j-i).c_str()));

        // set up for next loop
    i = j+1;
  }
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
      << "Multi-channel patch based segmentation tool" << std::endl
      << "Reference: Patch-based segmentation using expert priors: Application to hippocampus and ventricle segmentation" << std::endl
      << "\thttp://dx.doi.org/10.1016/j.neuroimage.2010.09.018" << std::endl << std::endl
      << "Usage: "<<name<<"  <input1> [input2] .. [inputN] " << std::endl
      << "\t--train <training_data> use this samples list" << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--mask <mask>" <<std::endl
      << "\t--patch <n>" <<std::endl
      << "\t--search <n>" <<std::endl
//       << "\t--dist <output> output minimal patch distance" <<std::endl
      << "\t--cls  <output> output closest patch id" <<std::endl
      << "\t--beta <f> scaling coefficient for minimal distance, default 0.125 to be compatible with Pierric's code"<<std::endl
      << "\t--scaling <f> use this map for scaling"<<std::endl
      << "\t--discrete <n> use for discrete labeling with n number of classes, including  0, default 2"<<std::endl
      << "\t--confidence <output> output confidence map"<<std::endl
      << "\t--adist <output> output mean search distance map"<<std::endl
      << "\t--grading <output> output image grading map"<<std::endl
      << "\t--iter <n> maximum number of iterations when preselection is used"<<std::endl
      << "\t--extract <n> extract label N from the label image first"<<std::endl
      << "\t--top <n> select only top n samples from each training set"<<std::endl
      << "\t--float store results as floats"<<  std::endl
      << "\t--short store results as short" <<  std::endl
      << "\t--byte  store results as  byte"  << std::endl
      << "\t--box use Box structuring element (default)" << std::endl
      << "\t--ball use ball structuring element "<< std::endl
      << "\t--prelabel <labels> prelabel data using this labels, segment only parts which are not labeled yet" << std::endl
      << "\t--prob <prefix> store per-label probabilities" << std::endl
      << "\t--output <output labels> store output " << std::endl
      << "\t--weights <w1,w2...>" << std::endl
      << "\t--groups <n> - modifies behaviour to allow one more column in the training sample (group id) and gurantees balanced sample after pre-selection"<<std::endl;
}

//! a helper function for image saving
template<class T> void save_image(typename T::Pointer image,const char *fname,
                                  const char* minc_datatype,
                                  itk::Object* metadata,
                                  const std::string& history)
{
  if(metadata) copy_metadata(image,metadata);
  if(!history.empty()) append_minc_history(image,history);
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

//! a helper function for multiple image loading
template <class T1,class T2> typename T1::Pointer load_images(const strings& files)
{
  typedef itk::ComposeImageFilter<T2>   ImageToVectorImageFilterType;
  typename ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
  for(size_t i=0;i<files.size();i++)
  {
    typename itk::ImageFileReader<T2>::Pointer reader = itk::ImageFileReader<T2>::New();
    reader->SetFileName(files[i].c_str());
    reader->Update();
    //if(!i) imageToVectorImageFilter->SetInput(reader->GetOutput());
    imageToVectorImageFilter->SetInput(i, reader->GetOutput());
  }
  imageToVectorImageFilter->Update();
  
  return imageToVectorImageFilter->GetOutput();
}


//! helper function that calculates SSD distances between all images in the train_images and in_img withing ROI
void ssd_distance(std::vector<ImagePairType>& train_images,
                  MultipleFeatureImageType::Pointer in_img,
                  itk::mask3d::Pointer roi_in,
                  bool verbose)
{
  if(verbose)
    std::cout<<"Calculating SSD.."<<std::endl;
  
  for(size_t i=0;i<train_images.size();i++)
  {
    //let's do it old fashioned way, for now....
    itk::ImageRegionIterator<MultipleFeatureImageType> 
        img_it(in_img, in_img->GetLargestPossibleRegion());
        
    itk::ImageRegionIterator<MultipleFeatureImageType> 
        sample_it(train_images[i].feature, in_img->GetLargestPossibleRegion());
        
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
          MultipleFeatureImageType::PixelType::RealValueType d=(img_it.Get()-sample_it.Get()).GetNorm();
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
        MultipleFeatureImageType::PixelType::RealValueType d=(img_it.Get()-sample_it.Get()).GetNorm();
        dist+=d*d;
        cnt++;
        ++img_it;
        ++sample_it;
      }
    }
    if(cnt>0) dist/=cnt;
    train_images[i].similarity=dist;
    
    if(verbose)
      std::cout<<"\tImage "<<i<<" SSD="<<dist<<std::endl;
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
  std::string history=minc_time_stamp(argc, argv); 
  int discrete=2;
  int max_iterations=50;
  int ball_structuring_element=0;
  std::string output_f;
  std::vector<double> weights;
  int groups=0;

  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"train", required_argument, 0, 't'},
    {"threshold", required_argument, 0, 'h'},
    {"mask", required_argument, 0, 'm'},
    {"patch", required_argument, 0, 'p'},
    {"search", required_argument, 0, 's'},
    {"beta", required_argument, 0, 'b'},
//     {"dist", required_argument, 0, 'd'},
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
    {"weights", required_argument, 0, 'w'},
    {"groups", required_argument, 0, 'g'},
    
    {"output", required_argument, 0, 'O'},
    
    {"box", no_argument, &ball_structuring_element, 0},
    {"ball", no_argument, &ball_structuring_element, 1},
    {0, 0, 0, 0}
  };

  int c;
  for ( ;; )
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long ( argc, argv, "t:r:m:T:D:S:C:G:I:E:P:o:O:w:", long_options, &option_index );

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
      case 'O':
        output_f=optarg;
        break;
      case 'w':
        stringtokd(weights,optarg,",");
        break;
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

  strings inputs_f;
  for(int i=optind;i<argc;i++)
    inputs_f.push_back(argv[i]);
  
  int channels=inputs_f.size();
  
  if(weights.empty())
  {
    for(int i=0;i<channels;i++)
      weights.push_back(1.0);
  } else if(weights.size()!=channels) {
    std::cerr << "Improper number of weights specifyed, expecting:"<<channels<<" got:"<<weights.size()<<std::endl;
    return 1;
  }
  

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
    MultipleFeatureImageType::Pointer in_img(MultipleFeatureImageType::New());
    SingleFeatureImageType::Pointer   out_conf,out_adist,out_grading;
    SingleFeatureImageType::Pointer   preload,preload_confidence;
    LabelImageType::Pointer   preload_discrete;
    
    
    in_img=load_images<MultipleFeatureImageType,SingleFeatureImageType>(inputs_f);

    itk::mask3d::Pointer roi_in(itk::mask3d::New());
    
    if(!mask_f.empty())
      roi_in=load_image<itk::mask3d>(mask_f.c_str());
    
    if(!preload_f.empty())
    {
      preload_discrete=load_image<LabelImageType>(preload_f.c_str());
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
    
    string_table tbl;
    
    read_table_n ( train_f.c_str(),tbl, channels+1 );
    if ( !table_verify ( tbl ) )
      return 1;
    
    if(verbose)
      std::cout<<"Loading volumes..."<<std::endl;
    
    std::vector<ImagePairType> train_images;
    std::vector<MultipleFeatureImageType::Pointer> trg_images;
    
    if(!groups)
    {
    
      for(size_t i=0;i<tbl.size();i++)
      {
          
        if(verbose) 
          std::cout<<tbl[i][0].c_str()<<"\t"<<std::flush;
        
        strings features_f;
        for(int j=0;j<channels;j++)
          features_f.push_back(tbl[i][j]);
        
        MultipleFeatureImageType::Pointer features=load_images<MultipleFeatureImageType,SingleFeatureImageType>(features_f);
        
        if(verbose) 
          std::cout<<tbl[i][1].c_str()<<"\t"<<std::flush;
        
        double grading=1;
        
        if(tbl[i].size()>(channels+1))
        {
          grading=atof(tbl[i][channels+1].c_str());
          if(verbose)
            std::cout<<"\t"<<grading;
        }
        if(verbose) 
          std::cout<<std::endl;

        SingleFeatureImageType::Pointer label=load_image<SingleFeatureImageType>(tbl[i][channels].c_str());
        train_images.push_back(ImagePairType(features, label, grading));
      }
      
      if(select_top)
      {
        // calculate SSD similarities
        ssd_distance(train_images ,in_img,roi_in,verbose);

        std::sort(train_images.begin(), train_images.end(), ssd_compare);

        if(train_images.size()>select_top)
          train_images.erase(train_images.begin()+select_top,train_images.end()); 

      }
    } else {
      std::set<int> group_set;
      for(size_t i=0;i<tbl.size();i++)
      {
        if(verbose) 
          std::cout<<tbl[i][0].c_str()<<"\t"<<std::flush;
          
        strings features_f;
        for(int j=0;j<channels;j++)
          features_f.push_back(tbl[i][j]);
        
        MultipleFeatureImageType::Pointer features=load_images<MultipleFeatureImageType,SingleFeatureImageType>(features_f);
        SingleFeatureImageType::Pointer label=load_image<SingleFeatureImageType>(tbl[i][channels].c_str());
        
        int    group=-1;
        double grading=1;
        
        if(tbl[i].size()>(channels+1))
        {
          grading=atof(tbl[i][channels+1].c_str());
          if(verbose)
            std::cout<<"\t"<<grading;
        }
        if(tbl[i].size()>(channels+2))
        {
          group=atoi(tbl[i][channels+1].c_str());
          if(verbose)
            std::cout<<"\t"<<group;
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
        train_images.push_back(ImagePairType(features,label,grading,0,group));
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
              <<"\tSmoothing parameter Beta="<<beta<<std::endl
              <<"\tUsing "<<(ball_structuring_element?"Ball":"Box")<<" structural element"<<std::endl;


    using SizeT=itk::Size<3>::SizeValueType;
    StructuringElementType::RadiusType   search_radius_ = {{static_cast<SizeT>(search_radius),static_cast<SizeT>(search_radius),static_cast<SizeT>(search_radius)}};
    StructuringElementType::RadiusType   patch_radius_  = {{static_cast<SizeT>(patch_radius),static_cast<SizeT>(patch_radius),static_cast<SizeT>(patch_radius)}};

    StructuringElementType  searchKernel=
      ball_structuring_element?StructuringElementType::Box(search_radius_):StructuringElementType::Ball(search_radius_);
      
    StructuringElementType  patchKernel=
      ball_structuring_element?StructuringElementType::Box(patch_radius_):StructuringElementType::Ball(patch_radius_);

    DiscreteLabelingType::Pointer discreteWeightedAverage=DiscreteLabelingType::New();

    discreteWeightedAverage->GetWeight().SetBeta(beta);
    
    discreteWeightedAverage->SetLabelCount(discrete);
    
    //DiscreteBaseFilterType::Pointer flt;
    CommandProgressUpdate::Pointer observer = CommandProgressUpdate::New();
    
    DiscreteSegmentationFilterNoPsType::Pointer filter(DiscreteSegmentationFilterNoPsType::New());
    
    filter->SetSegmentationLibrary(train_images);
    filter->SetSearchKernel(searchKernel);
    filter->SetPatchKernel(patchKernel);
    filter->SetProcess(discreteWeightedAverage);
    
    WeightedVectorPatchDistance::Pointer  weighted_distance=WeightedVectorPatchDistance::New();
    weighted_distance->set_weights(weights);
    filter->SetDistance(weighted_distance);
    
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
        save_image<SingleFeatureImageType>(filter->GetProbabilityMap(l),tmp,
                                          store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);            
      }
    }
    
    if(verbose)
      std::cout<<"done"<<std::endl;
    
    if(!output_f.empty())
      save_image<LabelImageType>(filter->GetOutput(),output_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);
    
    if(!output_conf_f.empty())
      save_image<SingleFeatureImageType>(out_conf,output_conf_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);
    
    if(!output_adist_f.empty())
      save_image<SingleFeatureImageType>(out_adist,output_adist_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);
    
    if(!output_grading_f.empty())
      save_image<SingleFeatureImageType>(out_grading,output_grading_f.c_str(),store_short?typeid(unsigned short).name():(store_byte?typeid(unsigned char).name():typeid(float).name()),in_img,history);

  } catch( itk::ExceptionObject & err ) {
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return 2;
  } 
  return 0;
}
