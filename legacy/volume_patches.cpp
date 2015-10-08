/* ----------------------------- MNI Header -----------------------------------
@NAME       :  volume_patches
@DESCRIPTION:  non-local patch based segmentation
@COPYRIGHT  :
              Copyright 2012 Vladimir Fonov, McConnell Brain Imaging Centre, 
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


#include "utils.h"
#include "minc_1_rw.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include "minc_1_simple.h"
#include "minc_1_simple_rw.h"
#include <getopt.h>
#include <math.h>
#include <unistd.h>

using namespace minc;
using namespace std;

const double epsilon=1e-5;

inline double lorentzian_p ( double x,double sigma )
{
  return log ( 1.0+0.5* ( x/sigma ) * ( x/sigma ) );
}

inline double lorentzian_psi ( double x,double sigma )
{
  return 2*x/ ( 2*sigma*sigma+x*x );
}
template<class T> double l2_patch_distance(
  const minc::simple_volume<T>& ima,
  const minc::simple_volume<T>& fima,
  int x,int y,int z,
  int xx,int yy,int zz,
  int f=1 )
{

  double distancetotal=0.0;
  int i,j,k;
  int ni1,nj1,ni2,nj2,nk1,nk2,kk;

  int s= ( 2*f+1 );

  for ( k=-f; k<=f; k++ ) {
    for ( j=-f; j<=f; j++ ) {
      for ( i=-f; i<=f; i++ ) {
        //if(!i && !j && !k) continue;
        ni1=x+i;
        nj1=y+j;
        nk1=z+k;

        ni2=x+xx+i;
        nj2=y+yy+j;
        nk2=z+zz+k;

        double d=ima.safe_get ( ni1,nj1,nk1 ) - fima.safe_get ( ni2,nj2,nk2 );

        distancetotal += d*d;
      }
    }
  }
  distancetotal/= ( s*s*s );
  return distancetotal;
}

void apply_patches ( const std::vector<volumes>& samples,
                     const minc_byte_volume& mask,
                     const minc_float_volume& input,
                     minc_float_volume& output,
                     minc_float_volume& dist,
                     minc_float_volume& cls,
                     double sigma,
                     int patch_radius,
                     int search_radius,
                     double beta,
                     bool progress=false )
{
  bool calc_labels=!output.empty();
  bool calc_dist  =!dist.empty();
  bool calc_cls   =!cls.empty();
  sigma*=sigma;

  int number_of_samples=samples.size()*(2*search_radius+1)*(2*search_radius+1)*(2*search_radius+1);

  std::vector<double> probe     ( number_of_samples,0 );
  std::vector<double> distance  ( number_of_samples,0 );
  std::vector<double> distance2 ( number_of_samples,0 );
  double min_distance;
  int    best_sample;

  std::cout<<"Sigma="<<sigma<<std::endl;
  std::cout<<"Beta="<<beta<<std::endl;

  if(calc_dist)
    std::cout<<"Calculating distances"<<std::endl;

  if(calc_labels)
    std::cout<<"Calculating labels"<<std::endl;

  if(calc_cls)
    std::cout<<"storing best sample"<<std::endl;

  std::vector<bool> progress_marks(20,false);

  for ( int z=0; z<input.dim ( 2 ); z++ ) {
    //if ( progress ) std::cout<<"."<<std::flush;
    if(progress) {
      if(!progress_marks[z*20/input.dim ( 2 )]) {
        std::cout<<z*100/input.dim(2)<<"% .. "<<std::flush;
        progress_marks[z*20/input.dim ( 2 )]=true;
      }
    }
    for ( int y=0; y<input.dim ( 1 ); y++ )
      for ( int x=0; x<input.dim ( 0 ); x++ ) {
        //gather statistics
        if ( !mask.get ( x,y,z ) ) {
          if(calc_cls)    cls.set ( x,y,z,-1 );
          if(calc_dist)   dist.set ( x,y,z,0.0 );
          if(calc_labels) output.set ( x,y,z,0.0 );
          continue;
        }

        int cnt=0;
        min_distance=1e10;
        for(int zz=-search_radius; zz<=search_radius; zz++) {
          for(int yy=-search_radius; yy<=search_radius; yy++) {
            for(int xx=-search_radius; xx<=search_radius; xx++) {
              for ( int i=0; i<samples.size(); i++ ) {
                distance[cnt]=l2_patch_distance<float> ( input,samples[i][0],x,y,z,xx,yy,zz,patch_radius );
                probe[cnt]   =samples[i][1].safe_get ( x+xx,y+yy,z+zz );

                if ( !cnt || distance[cnt] < min_distance ) {
                  min_distance=distance[cnt];
                  best_sample=cnt;
                }
                cnt++;
              }
            }
          }
        }
        /*
        distance2=distance;
        std::sort ( distance2.begin(),distance2.end() );

        double median_distance= ( number_of_samples&1 ) ?
                                ( sqrt ( distance2[number_of_samples/2] ) ) :
                                ( sqrt ( distance2[number_of_samples/2] ) +sqrt ( distance2[number_of_samples/2+1] ) ) /2;

        */
        if(calc_dist) dist.set ( x,y,z,min_distance );
        if(calc_cls)  cls.set ( x,y,z,best_sample+1 );

        //median_distance=min_distance;

        double out=0.0;
        double sum=0.0;
        /*
                median_distance*=median_distance;
                if ( median_distance<epsilon )
                  median_distance=epsilon;
        */
        if(calc_labels) {
          if( min_distance< epsilon )
            min_distance=epsilon;

          for ( int i=0; i<number_of_samples; i++ ) {

            double weight;
            if ( sigma>0.001 )
              weight =::exp ( -distance[i]/sigma );
            else
              weight =::exp ( -distance[i]/(beta*min_distance) );

            out+=probe[i]*weight;
            sum+=weight;
          }
          out/=sum;

          output.set ( x,y,z,out );
        }
      }
  }
  if(progress)
    std::cout<<std::endl;
}

void show_usage ( const char *name )
{
  std::cerr
      << "Usage: "<<name<<"  <input> [output labels]" << std::endl
      << "\t--train <training_data> use this samples list" << std::endl
      << "\t--train2 <training_data2> use this samples list" << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--mask <mask>" <<std::endl
      << "\t--patch <n>" <<std::endl
      << "\t--search <n>" <<std::endl
      << "\t--dist <output> output minimal patch distance" <<std::endl
      << "\t--cls  <output> output closest patch id" <<std::endl
      << "\t--beta <f> scaling coefficient for minimal distance"<<std::endl
      << "\t--scaling <f> use this map for scaling"<<std::endl;

}

int main ( int argc,char **argv )
{
  int clobber=0;
  int verbose=0;
  std::string train_f,train_f2;
  int selected;
  std::string mask_f;
  double sigma=0.0;
  int patch_radius=1;
  int search_radius=1;
  double beta=1.0;
  std::string scaling_f;
  std::string output_dist_f,output_cls_f;

  static struct option long_options[] = {
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"train", required_argument, 0, 't'},
    {"train2", required_argument, 0, 'T'},
    {"mask", required_argument, 0, 'm'},
    {"patch", required_argument, 0, 'p'},
    {"search", required_argument, 0, 's'},
    {"beta", required_argument, 0, 'b'},
    {"dist", required_argument, 0, 'd'},
    {"cls", required_argument, 0, 'c'},
    {"scaling", required_argument, 0, 'g'},
    {0, 0, 0, 0}
  };

  int c;
  for ( ;; ) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long ( argc, argv, "t:r:m:T:", long_options, &option_index );

    /* Detect the end of the options. */
    if ( c == -1 )
      break;

    switch ( c ) {
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
    case 'b':
      beta=atof(optarg);
      break;
    case 'd':
      output_dist_f=optarg;
      break;
    case 'c':
      output_cls_f=optarg;
      break;
    case 'g':
      scaling_f=optarg;
      break;
    case '?':
      /* getopt_long already printed an error message. */
    default:
      show_usage ( argv[0] );
      return 1;
    }
  }

  if ( ( argc - optind ) < 1 ) {
    show_usage ( argv[0] );
    return 1;
  }

  std::string input_f=argv[optind];
  std::string output_f;

  if(( argc - optind ) > 1)
    output_f=argv[optind+1];


  strings inputs_f ( 1 );

  inputs_f[0]=input_f;

  if (!clobber && !output_f.empty() && !access (output_f.c_str(), F_OK)) {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }

  if (!clobber && !output_dist_f.empty() && !access (output_dist_f.c_str(), F_OK)) {
    std::cerr << output_dist_f.c_str () << " Exists!" << std::endl;
    return 1;
  }

  if (!clobber && !output_cls_f.empty() && !access (output_cls_f.c_str(), F_OK)) {
    std::cerr << output_cls_f.c_str () << " Exists!" << std::endl;
    return 1;
  }

  try {
    string_table tbl;
    read_table_n ( train_f.c_str(),tbl,0 );

    //TODO: figure out if we need to adjust weights for the addditional samples
    if(!train_f2.empty()) {
      if(verbose)
        std::cout<<"Loading secondary training set!"<<std::endl;

      string_table tbl2;
      read_table_n ( train_f2.c_str(),tbl2,0 );

      for(int i=0; i<tbl2.size(); i++)
        tbl.push_back(tbl2[i]);
    }

    if ( !table_verify ( tbl ) )
      return 1;

    volumes inputs;
    std::vector<volumes> all_volumes;

    if(verbose)
      std::cout<<"Loading volumes..."<<std::flush;

    load_volumes ( inputs_f,inputs,0,verbose );


    load_all_volumes_0(tbl,all_volumes,tbl.size(),verbose);

    if(verbose)
      std::cout<<"done"<<std::endl;


    minc_byte_volume mask;
    minc_float_volume cls,dist,output;

    if ( !mask_f.empty() ) {
      minc_1_reader rdr;
      rdr.open ( mask_f.c_str() );
      load_simple_volume ( rdr,mask );
    } else {
      mask.resize ( inputs[0].size() );
      mask=1; //all is allowed
    }

    if(!output_cls_f.empty()) {
      cls.resize ( inputs[0].size() );
      cls=-1;
    }

    if(!output_dist_f.empty()) {
      dist.resize ( inputs[0].size() );
      dist=0.0;
    }

    if(!output_f.empty()) {
      output.resize ( inputs[0].size() );
      output=0.0;
    }


    //2. calculate local similarities
    if(verbose)
      std::cout<<"Patch based analysis, patch size="<<patch_radius<<" "
               <<"search radius="<<search_radius<<" "
               <<"samples="<<all_volumes.size()<<" ..."<<std::flush;

    apply_patches ( all_volumes,mask,inputs[0],output,dist,cls,sigma,patch_radius,search_radius,beta,verbose );

    if(verbose)
      std::cout<<"done"<<std::endl;


    minc_1_reader rdr;
    rdr.open ( input_f.c_str(),false,true );

    if(!output_f.empty()) {
      minc_1_writer wrt;
      wrt.open ( output_f.c_str(),rdr.info(),2,NC_FLOAT );
      save_simple_volume<float> ( wrt,output );
    }

    if(!output_dist_f.empty()) {
      minc_1_writer wrt;
      wrt.open ( output_dist_f.c_str(),rdr.info(),2,NC_FLOAT );
      save_simple_volume<float> ( wrt,dist );

    }

    if(!output_cls_f.empty()) {
      minc_1_writer wrt;
      wrt.open ( output_cls_f.c_str(),rdr.info(),2,NC_FLOAT );
      save_simple_volume<float> ( wrt,cls );
    }



  } catch ( const minc::generic_error & err ) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg() <<std::endl;
    return 1;
  }

  return 0;
}
