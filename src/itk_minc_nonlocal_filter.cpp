/* ----------------------------- MNI Header -----------------------------------
@NAME       :  itk_minc_nonlocal_filter
@DESCRIPTION:  non-local patch based filtering
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
#include "itkClassicalNonLocalFilter.h"

#include <itkFlatStructuringElement.h>
#include <itkSubtractImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkLogImageFilter.h>

#include "itkLocalMeanFilter.h"
#include "itkLocalSDFilter.h"
#include "itkMinimalDistanceNonLocalFilter.h"
#include "itkMedianDistanceNonLocalFilter.h"
#include "itkVariableNoiseNonLocalFilter.h"
#include "itkAdaptativeNonLocalFilter.h"
#include "itkNoiseMedianDistanceNonLocalFilter.h"
#include "itkMeanWeightNonLocalFilter.h"

#include <unistd.h>
#include <getopt.h>

#include <iostream>

//typedef itk::BinaryBallStructuringElement< itk::mask3d::PixelType, 3  >  FlatStructuringElementType;
typedef itk::FlatStructuringElement<3> FlatStructuringElementType;

typedef itk::image3d InputImageType;
typedef itk::image3d OutputImageType;


typedef itk::MinimalDistanceNonLocalFilter<
                            InputImageType,
                            OutputImageType,
                            FlatStructuringElementType,
                            FlatStructuringElementType,
                            itk::L2PatchDistance<InputImageType,FlatStructuringElementType> >  MinimalDistanceFilterType;

typedef itk::MedianDistanceNonLocalFilter<
                            InputImageType,
                            OutputImageType,
                            FlatStructuringElementType,
                            FlatStructuringElementType,
                            itk::L2PatchDistance<InputImageType,FlatStructuringElementType> >  MedianDistanceFilterType;

typedef itk::AdaptativeNonLocalFilter<
                            InputImageType,
                            OutputImageType,
                            FlatStructuringElementType,
                            FlatStructuringElementType,
                            itk::L2PatchDistance<InputImageType,FlatStructuringElementType>,
                            itk::InvExpWeight<double> >  AdaptativeFilterType;

typedef itk::NoiseMedianDistanceNonLocalFilter<
                            InputImageType,
                            OutputImageType,
                            FlatStructuringElementType,
                            FlatStructuringElementType,
                            itk::L2PatchDistance<InputImageType,FlatStructuringElementType>,
                            itk::InvExpWeight<double> >  NoiseMedianDistanceFilterType;
                            

typedef itk::MeanWeightNonLocalFilter<
                            InputImageType,
                            OutputImageType,
                            FlatStructuringElementType,
                            FlatStructuringElementType,
                            itk::L2PatchDistance<InputImageType,FlatStructuringElementType>,
                            itk::InvExpWeight<double> >  SimilarityFilterType;
                            
typedef itk::LocalSDFilter<InputImageType,
                            OutputImageType,
                            FlatStructuringElementType> SDFilterType;

typedef itk::LocalMeanFilter<InputImageType,
                            OutputImageType,
                            FlatStructuringElementType> MeanFilterType;
                            
                            
typedef itk::ImageToImageFilter<InputImageType,OutputImageType> ImageToImageFlt;

typedef itk::SubtractImageFilter<InputImageType,InputImageType,OutputImageType > SubtractImageFilterType;

#if (ITK_VERSION_MAJOR==3)
typedef itk::UnaryFunctorImageFilter<OutputImageType,OutputImageType,itk::Function::Log<OutputImageType::PixelType,OutputImageType::PixelType> > LogImageFlt;
#else
typedef itk::UnaryFunctorImageFilter<OutputImageType,OutputImageType,itk::Functor::Log<OutputImageType::PixelType,OutputImageType::PixelType> > LogImageFlt;
#endif

using namespace  std;

void show_usage (const char *name)
{
  std::cerr 
    << "Patch-based denoising tool" << std::endl
    << "Reference: Fast Non Local Means Denoising for 3D MR Images" << std::endl
    << "\thttp://dx.doi.org/10.1007/11866763_5" <<std::endl
    << "\tAdaptive non-local means denoising of MR images with spatially varying noise levels" << std::endl
    << "\thttp://dx.doi.org/10.1002/jmri.22003" << std::endl << std::endl
    << "Usage: " << name << " <input_image.mnc> <output_image.mnc> " << std::endl
    << "--noise <f> provide noise levels for filter that need it, use noise_estimate"<<std::endl
    << "--clobber clobber output files" << std::endl
    << "--verbose be more verbose" << std::endl
    << "--search <r> search radius in voxels, 1 - 3x3x3 search area, 2 - 5x5x5 etc default 2"<< std::endl
    << "--patch <r>  patch radius in voxels, 1 - 3x3x3 search area, 2 - 5x5x5 etc  default 1"<< std::endl
    << "--float save image in float format" << std::endl
    << "--short save image in short format" << std::endl
    << "--byte save image in byte format"   << std::endl
    << "--beta <f> beta parameter for Adaptative filter, default 1.0" << std::endl
    << "--mindist - calculate minimal patch distance " << std::endl
    << "--median - calculate median patch distance " << std::endl
    << "--similarity - local similarity metric" << std::endl
    << "--mweight  - calculate mean weight " << std::endl
    << "--mean - calculate local means " << std::endl
    << "--sd   - calculate local sd " << std::endl
    << "--flat - use flat structuring element " << std::endl
    << "--ball - use ball structuring element (default)" << std::endl
    << "--anlm - use adaptative non-local means filter " <<std::endl
    << "--roi <minc> - use ROI "<<std::endl
    << "--log - apply log transform to output image"<<std::endl
    << "--regularize <f> sigma for ANLM noise regularization kernel, default 0 (disabled)"<<std::endl
    << "--preselect - use preselection filter "<<std::endl
    << "--no_preselect - don't use preselection filter (default)"<<std::endl
    << "--th_m <f> - threshold for means preselection , default 0.95"<<std::endl
    << "--th_v <f> - threshold for means preselection , default 0.5"<<std::endl;
}

class CommandProgressUpdate : public itk::Command
{
  public:
    typedef CommandProgressUpdate   Self;
    typedef itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer; 
    itkNewMacro( Self );
    std::vector<bool>                _progress;
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
};


template<class TElement,
         class TDistance, 
         class TWeight,
         class TPreselect > 
typename itk::ClassicalNonLocalFilter<
                  InputImageType,
                  OutputImageType,
                  TElement,
                  TElement,
                  itk::L2PatchDistance<InputImageType,TElement>,
                  itk::InvExpWeight<double>,
                  TPreselect >::Pointer
    CreateClassicalNonLocalFilter(int search, int patch, double noise, 
                                  InputImageType::Pointer img_in,
                                  itk::mask3d::Pointer roi_in,
                                  bool calc_mweight, 
                                  bool flat_element )
{
    typedef itk::ClassicalNonLocalFilter<
                  InputImageType,
                  OutputImageType,
                  TElement,
                  TElement,
                  itk::L2PatchDistance<InputImageType,TElement>,
                  itk::InvExpWeight<double>,
                  TPreselect >  
          NonLocalFilterType;

    typename NonLocalFilterType::Pointer flt=NonLocalFilterType::New();
    itk::Size<3> srad; srad.Fill(search);
    itk::Size<3> prad; prad.Fill(patch);

    TElement  searchKernel;
    TElement  patchKernel;
    
    if(flat_element)
    {
      searchKernel =TElement::Box( srad);
      patchKernel = TElement::Box( prad );
    }  else {
      searchKernel=TElement::Ball( srad);
      patchKernel=  TElement::Ball( prad );
    }
    
    flt->SetSearchKernel( searchKernel );
    flt->SetPatchKernel( patchKernel );
    
    flt->Setsigma2( noise*noise );
    flt->SetInput(img_in);
    flt->SetOutputMeanWeight(calc_mweight);
    
    if(roi_in.IsNotNull())
      flt->SetRoiImage(roi_in);
    
    //flt->PrintSelf(std::cout,1);
    
    return flt;
}

template<class TElement,
         class TDistance, 
         class TWeight,
         class TPreselect > 
  ImageToImageFlt::Pointer
    CreateClassicalNonLocalFilterWithPS(int search, int patch, 
                                        double noise, InputImageType::Pointer img_in,
                                        itk::mask3d::Pointer roi_in, bool calc_mweight,
                                        bool flat_element,
                                        double th_m, double th_v )
{
  typedef itk::ClassicalNonLocalFilter<
                InputImageType,
                OutputImageType,
                TElement,
                TElement,
                itk::L2PatchDistance<InputImageType,TElement>,
                itk::InvExpWeight<double>,
                TPreselect >
        NonLocalFilterType;

  typename NonLocalFilterType::Pointer flt = CreateClassicalNonLocalFilter
        <TElement,TDistance,TWeight,TPreselect> (search,patch,noise,img_in,roi_in, calc_mweight, flat_element);
  
  typename TPreselect::Pointer _preselect=TPreselect::New();
  
    itk::Size<3> prad; prad.Fill(patch);

    TElement  patchKernel;
    
    if(flat_element)
    {
      patchKernel = TElement::Box( prad );
    }  else {
      patchKernel=  TElement::Ball( prad );
    }
  
  _preselect->SetPatchKernel(patchKernel);
  _preselect->SetImage(img_in);
  _preselect->SetThresholds(th_m, th_v);
  flt->SetPreselectionFilter(_preselect);
  //std::cout<<flt<<std::endl;
  
  return flt.GetPointer();
}


template<class TElement,
         class TDistance, 
         class TWeight,
         class TPreselect > 
typename itk::AdaptativeNonLocalFilter<
              InputImageType,
              OutputImageType,
              TElement,
              TElement,
              itk::L2PatchDistance<InputImageType,FlatStructuringElementType>,
              itk::InvExpWeight<double>,
              TPreselect >::Pointer
    CreateAdaptativeNonLocalFilter(int search, int patch, 
                                  double beta, double regularize,  
                                  InputImageType::Pointer img_in,
                                  itk::mask3d::Pointer roi_in,
                                  bool calc_mweight, 
                                  bool flat_element )
{
    typedef itk::AdaptativeNonLocalFilter<
              InputImageType,
              OutputImageType,
              TElement,
              TElement,
              itk::L2PatchDistance<InputImageType,FlatStructuringElementType>,
              itk::InvExpWeight<double>,
              TPreselect >  
          NonLocalFilterType;

    typename NonLocalFilterType::Pointer flt=NonLocalFilterType::New();
    itk::Size<3> srad; srad.Fill(search);
    itk::Size<3> prad; prad.Fill(patch);

    TElement  searchKernel;
    TElement  patchKernel;
    
    if(flat_element)
    {
      searchKernel =TElement::Box( srad);
      patchKernel = TElement::Box( prad );
    }  else {
      searchKernel=TElement::Ball( srad);
      patchKernel=  TElement::Ball( prad );
    }
    
    flt->SetSearchKernel( searchKernel );
    flt->SetPatchKernel( patchKernel );
    
    flt->SetBeta( beta );
    flt->SetRegularize( regularize );
    
    flt->SetInput(img_in);
    flt->SetOutputMeanWeight(calc_mweight);
    
    if(roi_in.IsNotNull())
      flt->SetRoiImage(roi_in);
    
    //flt->PrintSelf(std::cout,1);
    
    return flt;
}

template<class TElement,
         class TDistance, 
         class TWeight,
         class TPreselect > 
  ImageToImageFlt::Pointer
    CreateAdaptativeNonLocalFilterWithPS(int search, int patch, 
                                         double beta, double regularize, 
                                         InputImageType::Pointer img_in,
                                         itk::mask3d::Pointer roi_in, bool calc_mweight,
                                         bool flat_element,
                                         double th_m, double th_v )
{
  typedef itk::AdaptativeNonLocalFilter<
              InputImageType,
              OutputImageType,
              TElement,
              TElement,
              itk::L2PatchDistance<InputImageType,FlatStructuringElementType>,
              itk::InvExpWeight<double>,
              TPreselect >
        NonLocalFilterType;

  typename NonLocalFilterType::Pointer flt = CreateAdaptativeNonLocalFilter
      <TElement,TDistance,TWeight,TPreselect> (search, patch, beta, 
                            regularize, img_in, roi_in, calc_mweight, flat_element);

  typename TPreselect::Pointer _preselect=TPreselect::New();

  itk::Size<3> prad; prad.Fill(patch);

  TElement  patchKernel;

  if(flat_element)
  {
    patchKernel = TElement::Box( prad );
  }  else {
    patchKernel=  TElement::Ball( prad );
  }
  
  _preselect->SetPatchKernel(patchKernel);
  _preselect->SetImage(img_in);
  _preselect->SetThresholds(th_m, th_v);
  flt->SetPreselectionFilter(_preselect);
  //std::cout<<flt<<std::endl;
  
  return flt.GetPointer();
}
                            
                            
                            
                            
                            
                            


int main (int argc, char **argv)
{
  int clobber=0;
  int verbose=0;
  int store_float=0;
  int store_short=0;
  int store_byte=0;
  int c;
  int search=2;
  int patch=1;
  int flat_element=1;
  int calc_mean=0;
  int calc_sd=0;
  int calc_mdist=0;
  int calc_anlm=0;
  int calc_median=0;
  int calc_log=0;
  int calc_similarity=0;
  int calc_mweight=0;
  
  int use_preselect=1;
  
  double th_m=0.95;
  double th_v=0.5;
  
  double noise=0.0;
  double beta=1.0;
  double anlm_regularize=0.0;
  
  int iterations=0;
  std::string history=itk::minc_time_stamp(argc, argv); 
  std::string roi_f;

  static struct option long_options[] = { 
    {"clobber", no_argument, &clobber,           1},
    {"verbose", no_argument, &verbose,           1},
    {"preselect", no_argument, &use_preselect,   1},
    {"no_preselect", no_argument, &use_preselect, 0},
    {"mean",    no_argument, &calc_mean,   1},
    {"sd",      no_argument, &calc_sd,     1},
    {"flat",    no_argument, &flat_element,1},
    {"ball",    no_argument, &flat_element,0},
    {"float",   no_argument, &store_float, 1},
    {"short",   no_argument, &store_short, 1},
    {"byte",    no_argument, &store_byte,  1},
    {"mindist",   no_argument,&calc_mdist,  1},
    {"mweight",   no_argument,&calc_mweight,  1},
    {"anlm",    no_argument, &calc_anlm,  1},
    {"median",  no_argument, &calc_median,  1},
    {"similarity",  no_argument, &calc_similarity,  1},
    {"log",     no_argument, &calc_log,  1},
    {"search",  required_argument, 0,    's'},
    {"patch",   required_argument, 0,    'p'},
    {"roi",     required_argument, 0,    'r'},
    {"noise",   required_argument, 0,    'n'},
    {"sigma",   required_argument, 0,    'n'},
    {"beta",    required_argument, 0,    'b'},
    {"regularize",required_argument, 0,    'e'},
    {"th_m",required_argument, 0,    'M'},
    {"th_v",required_argument, 0,    'V'},
    /*{"weights", required_argument, 0, 'w'},
    {"sweights", required_argument, 0, 'e'},*/
    {0, 0, 0, 0}};

  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "s:p:r:n:b:e:M:V:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
    case 0:
      break;
    case 's':
      search = atoi(optarg);
      break;
    case 'p':
      patch = atoi(optarg);
      break;
    case 'r':
      roi_f = optarg;
      break;
    case 'n':
      noise=atof(optarg);
      break;
    case 'b':
      beta=atof(optarg);
      break;
    case 'e':
      anlm_regularize=atof(optarg);
      break;
    case 'V':
      th_v=atof(optarg);
      break;
    case 'M':
      th_m=atof(optarg);
      break;
/*		case 'w':
      weights_f = optarg;
      break;
    case 'e':
      sweights_f = optarg;
      break;*/
    case '?':
    default:
      show_usage (argv[0]);
      return 1;
    }
  }

  if((argc - optind)<2)
  {
    show_usage(argv[0]);
    return 1;
  }
  
  std::string input(argv[optind]),output(argv[optind+1]);
  
  if (!clobber && !access (output.c_str (), F_OK))
  {
    cerr << output.c_str () << " Exists!" << endl;
    return 1;
  }
  
  try
  {
    itk::image3d::Pointer img_in;
    itk::image3d::Pointer out_img;
    
    img_in=itk::load_volume<itk::image3d>(input.c_str());
    
    itk::mask3d::Pointer roi_in;
    if(!roi_f.empty())
      roi_in=itk::load_volume<itk::mask3d>(roi_f.c_str());
    
    CommandProgressUpdate::Pointer observer = CommandProgressUpdate::New();
    ImageToImageFlt::Pointer filter;
    ImageToImageFlt::Pointer filter_1,filter_2,filter_3;
    ImageToImageFlt::Pointer log_filter;
    
    itk::Size<3> srad; srad.Fill(search);
    itk::Size<3> prad; prad.Fill(patch);

    FlatStructuringElementType  searchKernel=FlatStructuringElementType::Ball( srad );
    FlatStructuringElementType  patchKernel =FlatStructuringElementType::Ball( prad ) ;
    
    if( calc_sd )
    { 
      SDFilterType::Pointer flt(SDFilterType::New()); 
      FlatStructuringElementType  patchKernel;
      patchKernel.Ball( prad );
      
      flt->SetKernel( patchKernel );
      flt->SetInput(img_in);
      filter=flt;
    } else if( calc_mean ) {
      MeanFilterType::Pointer flt(MeanFilterType::New()); 
      FlatStructuringElementType  patchKernel;
      patchKernel.Ball( prad );
      
      flt->SetKernel( patchKernel );
      flt->SetInput(img_in);
      filter=flt;
      
    } else if( calc_mdist ) {
      MeanFilterType::Pointer fltm(MeanFilterType::New()); 
      MinimalDistanceFilterType::Pointer flt(MinimalDistanceFilterType::New());
      
      fltm->SetKernel( patchKernel );
      fltm->SetInput(img_in);
      //fltm->Update();
      SubtractImageFilterType::Pointer sub(SubtractImageFilterType::New());
      sub->SetInput1(img_in);
      sub->SetInput2(fltm->GetOutput());
      
      filter_1=fltm;
      filter_2=sub;
      
      flt->SetSearchKernel( searchKernel );
      flt->SetPatchKernel( patchKernel );
      flt->SetInput(sub->GetOutput());
      
      if(!roi_f.empty())
        flt->SetRoiImage(roi_in);
      
      filter=flt;
      
    } else if(calc_median){
      MedianDistanceFilterType::Pointer flt(MedianDistanceFilterType::New());
            
      flt->SetSearchKernel( searchKernel );
      flt->SetPatchKernel( patchKernel );
      
      flt->SetInput(img_in);
      
      if(!roi_f.empty())
        flt->SetRoiImage(roi_in);
      
      filter=flt;
    } else if(calc_similarity) {
      SimilarityFilterType::Pointer flt(SimilarityFilterType::New());
            
      flt->SetSearchKernel( searchKernel );
      flt->SetPatchKernel( patchKernel );
      flt->Setsigma2( noise*noise );

      flt->SetInput(img_in);
      
      if(!roi_f.empty())
        flt->SetRoiImage(roi_in);
      
      filter=flt;
    }	else if(calc_anlm) { // running adaptative filter
    
    
        if( use_preselect )
        {
          filter=CreateAdaptativeNonLocalFilterWithPS<
                  FlatStructuringElementType,
                  itk::L2PatchDistance<InputImageType,FlatStructuringElementType>,
                  itk::InvExpWeight<double>, itk::MeanAndSdPreselectionFilter<InputImageType,FlatStructuringElementType> >
              (search, patch, beta, anlm_regularize, img_in, roi_in, calc_mweight, flat_element, th_m, th_v);
        } else {
        
          filter=CreateAdaptativeNonLocalFilter<
                  FlatStructuringElementType,
                  itk::L2PatchDistance<InputImageType,FlatStructuringElementType>,
                  itk::InvExpWeight<double>, itk::NOOPPreselection<3> >
              (search, patch, beta, anlm_regularize, img_in, roi_in, calc_mweight, flat_element);
        }
      
    } else  { // running classic non-local patch filter
        if( use_preselect )
        {
          filter=CreateClassicalNonLocalFilterWithPS<
                  FlatStructuringElementType,
                  itk::L2PatchDistance<InputImageType,FlatStructuringElementType>,
                  itk::InvExpWeight<double>, itk::MeanAndSdPreselectionFilter<InputImageType,FlatStructuringElementType> >
              (search, patch, noise, img_in, roi_in, calc_mweight, flat_element, th_m, th_v);
        } else {
        
          filter=CreateClassicalNonLocalFilter<
                  FlatStructuringElementType,
                  itk::L2PatchDistance<InputImageType,FlatStructuringElementType>,
                  itk::InvExpWeight<double>, itk::NOOPPreselection<3> >
              (search, patch, noise, img_in, roi_in, calc_mweight, flat_element);
        }
    }
    if(filter_1.IsNotNull())
      filter_1->AddObserver(itk::ProgressEvent(),observer);
    
    if(filter_2.IsNotNull())
      filter_2->AddObserver(itk::ProgressEvent(),observer);
    
    if(filter_3.IsNotNull())
      filter_3->AddObserver(itk::ProgressEvent(),observer);

    filter->AddObserver(itk::ProgressEvent(),observer);

    if(calc_log)
    {
      log_filter=LogImageFlt::New();
      log_filter->AddObserver(itk::ProgressEvent(),observer);
      log_filter->SetInput(filter->GetOutput());
      
      log_filter->Update();
      out_img=log_filter->GetOutput();
    } else {
    
      filter->Update();
      out_img=filter->GetOutput();
    }
    
    itk::copy_metadata(out_img,img_in);
    itk::append_minc_history(out_img,history.c_str());

    if(store_float)
      itk::set_minc_storage_type(out_img,typeid(float).name());
    else if(store_short) 
      itk::set_minc_storage_type(out_img,typeid(unsigned short).name());
    else if(store_byte) 
      itk::set_minc_storage_type(out_img,typeid(unsigned char).name());

    itk::save_volume<itk::image3d>(output.c_str(),  out_img);
    
  } catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return 2;
  } 
  return 0;

}

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2
