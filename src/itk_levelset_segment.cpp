#include "itkHelpers.h"

#include <itkThresholdSegmentationLevelSetFunction.h>
#include <itkSegmentationLevelSetImageFilter.h>
#include <itkThresholdSegmentationLevelSetImageFilter.h>

#include <unistd.h>
#include <iostream>
#include <getopt.h>


typedef itk::ThresholdSegmentationLevelSetFunction< minc::image3d, minc::image3d> ThresholdLevelSetFunction;
typedef itk::ThresholdSegmentationLevelSetImageFilter< minc::image3d, minc::image3d, minc::image3d::PixelType> MyLevelSetFilter;
using namespace  std;
using namespace  minc;

void show_usage (const char *name)
{
  std::cerr 
    << "Minc wrapper around itk::ThresholdSegmentationLevelSetFunction and SegmentationLevelSetImageFilter"<<std::endl
    << "Usage: " << name << " <input_image.mnc> <input_segmentation.mnc> <output_segmentation.mnc> " << std::endl
    << "--clobber clobber output files" << std::endl
    << "--iter <n> - maximum number of iterations, default 0 (infinite)"<<std::endl
    << "--min <min>, default auto" << std::endl
    << "--max <max>, default auto" << std::endl
    << "--float save image in float format" << std::endl
    << "--short save image in short format" << std::endl
    << "--byte save image in byte format"   << std::endl;
}

class CommandProgressUpdate : public itk::Command
{
  public:
    typedef CommandProgressUpdate   Self;
    typedef itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer; 
    itkNewMacro( Self );
    std::vector<bool> 							_progress;
  protected:
    CommandProgressUpdate():
      _progress(101,false)
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
      int prg=filter->GetProgress()*100;
      if(!_progress[prg])
      {
        std::cout << prg << "% " <<std::flush;
        _progress[prg]=true;
      } else {
        std::cout << "."<<std::flush;
      }
    }
};

int main (int argc, char **argv)
{
  
  int clobber=0;
  int store_float=0;
  int store_short=0;
  int store_byte=0;
  int c;
  bool set_min_max=false;
  double feature_min=0.0,feature_max=0.0;
  
  int iterations=0;
#ifdef HAVE_MINC4ITK  
  char *_history = time_stamp(argc, argv); 
  std::string history=_history;
  free(_history);
#else
  std::string history=minc_timestamp(argc,argv);
#endif  

  static struct option long_options[] = { 
    {"clobber", no_argument, &clobber, 1},
    {"float", no_argument, &store_float, 1},
    {"short", no_argument, &store_short, 1},
    {"byte", no_argument, &store_byte, 1},
    {"iter", required_argument, 0, 'i'},
    {"min", required_argument, 0, 'm'},
    {"max", required_argument, 0, 'M'},
    {0, 0, 0, 0}};

  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "i:m:M:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
    case 0:
      break;
    case 'M':
      feature_max = atof(optarg);
      set_min_max=true;
      break;
    case 'm':
      feature_min = atof(optarg);
      set_min_max=true;
      break;
    case 'i':
      iterations = atoi(optarg);
      break;
    case '?':
    default:
      show_usage (argv[0]);
      return 1;
    }
  }

  if((argc - optind)<3)
  {
    show_usage(argv[0]);
    return 1;
  }
  
  std::string input_feature(argv[optind]),input_set(argv[optind+1]),output(argv[optind+2]);
  
  if (!clobber && !access (output.c_str (), F_OK))
  {
    cerr << output.c_str () << " Exists!" << endl;
    return 1;
  }
  
  try
  {
    image3d::Pointer img_in(image3d::New());
    image3d::Pointer mask_in(image3d::New());
    load_minc (input_feature.c_str(), img_in);
    load_minc (input_set.c_str(), mask_in);
    
    MyLevelSetFilter::Pointer fltSegment(MyLevelSetFilter::New());
    ThresholdLevelSetFunction::Pointer fnSegment(ThresholdLevelSetFunction::New());
    itk::Size<3> rad; rad.Fill(1);
    
    /*fnSegment->SetCurvatureWeight(1.0);
    fnSegment->SetAdvectionWeight(1.0);
    fnSegment->SetPropagationWeight(1.0);
    
    fnSegment->Initialize(rad);
    
    fnSegment->SetUpperThreshold(feature_max);
    fnSegment->SetLowerThreshold(feature_min);
    fnSegment->SetFeatureImage(img_in);
    fnSegment->CalculateSpeedImage();*/
    
  // Set the inputs to the segmentation filter
    fltSegment->SetUpperThreshold(feature_max);
    fltSegment->SetLowerThreshold(feature_min);
    //fltSegment->SetSegmentationFunction(fnSegment);
    
    fltSegment->SetIsoSurfaceValue(0);//?
    fltSegment->SetInput(mask_in); 
    fltSegment->SetFeatureImage(img_in);
    
    fltSegment->SetNumberOfLayers(3);
    //fltSegment->SetIsoSurfaceValue(0.0);
    
    if(iterations>0)
      fltSegment->SetNumberOfIterations(iterations);
    
    CommandProgressUpdate::Pointer observer = CommandProgressUpdate::New();

    fltSegment->AddObserver(itk::ProgressEvent(),observer);
    fltSegment->Update();
    
    typedef itk::ShiftScaleImageFilter<MyLevelSetFilter::OutputImageType, image3d> DummyFilter;
    DummyFilter::Pointer fltDummy = DummyFilter::New();
    fltDummy->SetInput(fltSegment->GetOutput());
    fltDummy->SetScale(1.0);
    fltDummy->SetShift(0.0);
    fltDummy->Update();
    
    
    image3d::Pointer out_img(fltDummy->GetOutput());
    
    copy_metadata(out_img,img_in);
    append_history(out_img,history);
#ifdef HAVE_MINC4ITK    
    if(store_float)
    {
      minc::set_minc_storage_type(out_img,NC_FLOAT,true);
    } else if(store_short) {
      minc::set_minc_storage_type(out_img,NC_SHORT,true);
    } else if(store_byte) {
      minc::set_minc_storage_type(out_img,NC_BYTE,false);
    }
#else
  if(store_float)
  {
    minc::set_minc_storage_type(out_img,typeid(float).name());
  } else if(store_short) {
    minc::set_minc_storage_type(out_img,typeid(unsigned short).name());
  } else if(store_byte) {
    minc::set_minc_storage_type(out_img,typeid(unsigned char).name());
  }
#endif
    
    save_minc<image3d>(output.c_str(),  out_img);
    
  } catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    cerr << err.msg()<<endl;
    return 1;
  }
  return 0;
  
}
// kate: space-indent on; hl c++;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 2 
