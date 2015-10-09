#! /usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;

my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $mydir=dirname($0);
my $model_dir;
my $nuc;
my $cleanup;
my $subject;
my $qc=1;
my $manual_xfm;
my $nuc;
my $mri_3t;
my $exclude;
my $select=50;
my $preprocess;
my $ad_lib;
my $nc_lib;
my $variant='';
my $use_mmse=0;

GetOptions (
          'verbose'        => \$verbose,
          'clobber'        => \$clobber,
          'nuc'            => \$nuc,
          'model-dir=s'    => \$model_dir,
          'subject=s'      => \$subject,
          '3t'             => \$mri_3t,
          'exclude=s'      => \$exclude,
          'select=n'       => \$select,
          'preprocess'     => \$preprocess,
          'ad_lib=s'       => \$ad_lib,
          'nc_lib=s'       => \$nc_lib,
          'variant=s'      => \$variant,
          'manual=s'       => \$manual_xfm,
          'use_mmse'       => \$use_mmse,
          );

my $HELP=<<HELP ;
Reference: Pierrick Coupé, Simon F Eskildsen, José V Manjón, Vladimir S Fonov, D Louis Collins, Alzheimer's disease Neuroimaging Initiative 
          "Simultaneous segmentation and grading of anatomical structures for patient's classification: application to Alzheimer's disease" 
          http://dx.doi.org/10.1016/j.neuroimage.2011.10.080

Usage $me input_t1.mnc output_prefix
          --model-dir <anatomical model directory> - required!
         [
          --subject <id> - subject id
          --nuc - perform preprocessing - Denoise + N3 + intensity normalization
          --3t - subject was scanned on a 3T scanner
          --exclude <s> - exclude this subject from the library
          --select <n> select top N subjects from each library
          --preprocess - perform only preprocessing (no segmentation)
          --ad_lib <ad library.lst>
          --nc_lib <nc library.lst>
          --variant <v>
          --manual <xfm> provide manual initialization for linear transformation
          --use_mmse - use MMSE scores instead of grading
         ]
HELP


die $HELP if $#ARGV<1 || !$model_dir;


my $model     ="$model_dir/ad_model/model_t1w.mnc";
my $model_mask="$model_dir/ad_model/model_t1w_mask.mnc";
my $model_bbox="$model_dir/ad_model/model_t1w_bbox.mnc";
my $model_outline="$model_dir/ad_model/model_t1w_outline.mnc";
my $model_snipe_mask="$model_dir/ad_model/model_t1w_mask_snipe.mnc";
my $model_snipe_mask_bbox="$model_dir/ad_model/model_t1w_mask_snipe_bbox.mnc";
my $model_hcec_mask="$model_dir/ad_model/model_hcec_mask.mnc";
my $model_hcec_mask_bbox="$model_dir/ad_model/model_hcec_mask_bbox.mnc";
my $model_outline_lle="$model_dir/ad_model/model_t1w_bbox_lle.mnc";

my $NC_library="$model_dir/NCQCLibrary3_mnc.lst";
my $AD_library="$model_dir/ADQCLibrary3_mnc.lst";

my $mmse_scale=1.0;
my $mmse_bias=0.0;

if($use_mmse)
{
  $NC_library="$model_dir/NCQCLibrary3_mnc_mmse.lst";
  $AD_library="$model_dir/ADQCLibrary3_mnc_mmse.lst";
  # calculated based on distribution of selected subjects
  $mmse_scale=3.288719;
  $mmse_bias=26.425;
}

$NC_library=$nc_lib if $nc_lib;
$AD_library=$ad_lib if $ad_lib;

my $NC_library_ex="${NC_library}_${subject}.lst";
my $AD_library_ex="${AD_library}_${subject}.lst";

$variant="_${variant}" if $variant;

check_presence($model,$model_mask,$model_bbox,$model_outline,$model_snipe_mask,$NC_library,$AD_library,$model_hcec_mask);

my ($in,$prefix)=@ARGV;
my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );
do_cmd('mkdir','-p',$prefix);

unless($subject) #figure out subject id
{
  $subject=basename($in,'.gz');
  $subject =~ s/\.gz$//;
  $subject =~ s/\.mnc$//;
}

my $compress=$ENV{MINC_COMPRESS};
#$ENV{ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}=1;

if($nuc) # run N3 first
{
  if(! -e "${prefix}/nuc_${subject}.mnc")
  {
    N3_correct($in,"${prefix}/nuc_${subject}.mnc",$model,$model_mask,$mri_3t);
  }
  $in="${prefix}/nuc_${subject}.mnc";
}


#linear registration
if(! -e "${prefix}/stx_${subject}.xfm")
{  
  if( $manual_xfm )
  {
    do_cmd('cp',$manual_xfm,"${prefix}/manual_${subject}.xfm");
  }
  
  if( -e "${prefix}/manual_${subject}.xfm")
  {
    do_cmd('bestlinreg_s',
          $in,$model,'-nmi',
          "${prefix}/stx_${subject}.xfm",
           '-init_xfm',"${prefix}/manual_${subject}.xfm");
  } else {
    do_cmd('bestlinreg_s',
          $in,$model,'-nmi',
          "${prefix}/stx_${subject}.xfm");
  }
}

if(! -e "${prefix}/stx_${subject}.mnc")
{
  do_cmd('itk_resample',
         $in,'--like',$model,
         '--transform',"${prefix}/stx_${subject}.xfm",
         "${prefix}/stx_${subject}.mnc"
        );
}

if($qc && ! -e "${prefix}/qc_stx_${subject}.jpg")
{
  do_cmd('minc_qc.pl',"${prefix}/stx_${subject}.mnc",
         '--mask',$model_outline,"${prefix}/qc_stx_${subject}.jpg",'--image-range',0,100);
}

# local linear registration
if(! -e "${prefix}/stx_${subject}_bbox.xfm")
{
  do_cmd('bestlinreg_s','-close',"${prefix}/stx_${subject}.mnc",$model,"$tmpdir/stx_${subject}_bbox.xfm",
         '-lsq9','-target_mask',$model_snipe_mask);

  do_cmd('xfmconcat',"${prefix}/stx_${subject}.xfm","$tmpdir/stx_${subject}_bbox.xfm","${prefix}/stx_${subject}_bbox.xfm");
}

if(! -e "${prefix}/stx_${subject}_bbox.mnc")
{
  do_cmd('itk_resample',
         $in,'--like',$model_bbox,
         '--transform',"${prefix}/stx_${subject}_bbox.xfm",
         "${prefix}/stx_${subject}_bbox.mnc"
        );
}

# second intensity normalization
if(! -e "${prefix}/stx_${subject}_bbox_snipe.mnc")
{
  # TODO: figure out optimal parameters
  do_cmd('minc_nuyl',"${prefix}/stx_${subject}_bbox.mnc",$model_bbox,"${prefix}/stx_${subject}_bbox_snipe.mnc",
         '--source-mask',$model_snipe_mask_bbox,
         '--target-mask',$model_snipe_mask_bbox);
}

if(! -e "${prefix}/qc_bbox_${subject}.jpg")
{
  do_cmd('minc_qc.pl',"${prefix}/stx_${subject}_bbox_snipe.mnc",
         '--mask',$model_outline_lle,"${prefix}/qc_bbox_${subject}.jpg",
         '--image-range',0,100);
}

die "Only preprocessing is needed!\n" if $preprocess;

# finally run SNIPE
my @labels=(2,4,21,19); #VF: fixed EC left-right swap
my @structures=('HC_right','HC_left','EC_right','EC_left');

my @patch= (3,3,2,2);
my @search=(4,4,4,4);

if($exclude)
{
  do_cmd("fgrep -v $exclude $NC_library>$NC_library_ex");
  do_cmd("fgrep -v $exclude $AD_library>$AD_library_ex");
  $NC_library=$NC_library_ex;
  $AD_library=$AD_library_ex;
}

my $i;

for($i=0;$i<=$#labels;$i+=1)
{
  my $l=$labels[$i];
  my $p=$patch[$i];
  my $s=$search[$i];
  my $mask="$model_dir/mask_${l}.mnc";


  do_cmd('itk_patch_morphology',
          "${prefix}/stx_${subject}_bbox_snipe.mnc",
          "$prefix/stx_${subject}_seg_label_${l}.mnc",
         '--grad',"$prefix/stx_${subject}_grad_label_${l}${variant}.mnc",
         '--extract',$l,
         '--patch',$p,'--search',$s,
         '--mask',$mask,
         '--threshold',0.97,
         '--top',$select,
         '--train',$NC_library,
         '--train2',$AD_library,
         '--verbose','--iter','5') 
    unless -e "$prefix/stx_${subject}_seg_label_${l}${variant}.mnc";
}

unless( -e "$prefix/stx_${subject}_HCEC${variant}.mnc")
{
 # need to threshold at 50% first
 for($i=0;$i<=$#labels;$i+=1)
 {
   my $l=$labels[$i];
   do_cmd('minccalc',
          '-express','A[0]>0.5?A[0]:0',
          "$prefix/stx_${subject}_seg_label_${l}${variant}.mnc",
          "$tmpdir/stx_${subject}_seg_label_${l}${variant}.mnc");
 }

 do_cmd('crispify', '-volume', "$prefix/stx_${subject}_HCEC${variant}.mnc",
  "$tmpdir/stx_${subject}_seg_label_4.mnc", 4,
  "$tmpdir/stx_${subject}_seg_label_2.mnc", 2,
  "$tmpdir/stx_${subject}_seg_label_21.mnc", 21,
  "$tmpdir/stx_${subject}_seg_label_19.mnc", 19);
}

if($qc && ! -e "${prefix}/qc_HCEC_${subject}.jpg")
{
  do_cmd('minc_qc.pl',"${prefix}/stx_${subject}_bbox_snipe.mnc",
         '--mask',"$prefix/stx_${subject}_HCEC.mnc","${prefix}/qc_HCEC_${subject}.jpg",
         '--image-range',0,100,'--labels-mask');
}


unless(-e "${prefix}/grading_${subject}.txt")
{
  open OUT,">${prefix}/grading_${subject}.txt" or die "Can't open ${prefix}/grading_${subject}.txt for writing!\n";

  print OUT "${subject}";

  for($i=0;$i<=$#labels;$i+=1)
  {
    my $l=$labels[$i];
    my $s=$structures[$i];
    my $grad=`mincstats -q -mean $prefix/stx_${subject}_grad_label_${l}.mnc -mask $prefix/stx_${subject}_seg_label_${l}.mnc -mask_binvalue 1`;
    chomp $grad;
    $grad=$grad*1.0;
    
    # converting to real MMSE value if needed
    $grad=$grad*$mmse_scale+$mmse_bias;
    
    print OUT ",$grad";
  }

  for($i=0;$i<=$#labels;$i+=1)
  {
    my $l=$labels[$i];
    my $s=$structures[$i];
    my $vol=`mincstats -q -volume $prefix/stx_${subject}_grad_label_${l}.mnc -mask $prefix/stx_${subject}_seg_label_${l}.mnc -mask_binvalue 1`;
    chomp $vol;
    print OUT ",$vol";
  }
  
  # calculate SNIPE volume scale
  my $_scaling=`xfm2param ${prefix}/stx_${subject}_bbox.xfm |fgrep scale|tr -s ' ' ','`;
  chomp($_scaling);
  my @scaling=split(',',$_scaling);
  my $scale=$scaling[1]*$scaling[2]*$scaling[3];
  print OUT ",$scale";
  print OUT "\n";
  close OUT;
}

if($exclude)
{
  unlink($NC_library_ex);
  unlink($AD_library_ex);
}


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
}
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}

sub check_presence {
  my $i;
  foreach $i(@_)
  {
    die("$i does not exists!\n") unless -e $i;
  }
}

sub N3_correct {
  my ($in,$out,$model,$model_mask,$mri_3t)=@_;
  
  do_cmd('cp', $in, "$tmpdir/in.mnc");
  #fix some broken minc files
  
  my $zspacing=`mincinfo -attvalue zspace:spacing $tmpdir/in.mnc`;
  chomp($zspacing);
  
  if($zspacing =~ /irregular/)
  {
    do_cmd('minc_modify_header','-sinsert',
           'zspace:spacing=regular__',"$tmpdir/in.mnc");
  }

  do_cmd('minc_anlm',"$tmpdir/in.mnc","$tmpdir/denoise.mnc");

  my @args=( "nu_correct", "-clobber", 
            "-iter", 100, 
            "-stop", 0.0001, 
            "-fwhm", 0.1,
            "$tmpdir/denoise.mnc",  "$tmpdir/nuc.mnc",
            '-clobber');

  push @args,'-distance',50 if $mri_3t;
  do_cmd(@args);

  do_cmd('volume_pol', '--order', 1, 
        '--min', 0, '--max', 100, 
        '--noclamp',
        "$tmpdir/nuc.mnc", $model,'--expfile', "$tmpdir/stats", '--clobber');
        
  do_cmd('minccalc', "$tmpdir/nuc.mnc", $out, 
        '-expfile', "$tmpdir/stats", '-clobber','-short');
}

# kate: space-indent on; indent-width 2; replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl perl
