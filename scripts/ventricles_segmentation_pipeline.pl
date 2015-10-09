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
my $nuyl;
my $whole_head;
my $subject;
my $nlmask;
my $qc=0;
my $manual_xfm;
my $search=1;
my $reference; #for longitudinal
my $nonbrainconstraint;
my $corr=1.0;
my $translation=0.0;

GetOptions (
          'verbose'        => \$verbose,
          'qc'             => \$qc,
          'clobber'        => \$clobber,
          'model-dir=s'    => \$model_dir,
          'nuc'            => \$nuc,
          'cleanup'        => \$cleanup,
          'nuyl'           => \$nuyl,
          'whole-head'     => \$whole_head,
          'subject=s'      => \$subject,
          'nlmask'         => \$nlmask,
          'manual-xfm=s'   => \$manual_xfm,
          'search=n'       => \$search,
          'reference=s'    => \$reference,
          'nonbrainconstraint' => \$nonbrainconstraint
          );

my $HELP=<<HELP ;
Usage $me input_t1.mnc output_prefix  
          --model-dir <anatomical model directory> - required!
         [--nuc  run N3 before doing other steps
          --cleanup remove intermediate files 
          --nuyl use Nuyl's normalization instead of linear
          --while-head use whole head nonlinear registration - will require BET !
          --subject <subject id>
          --nlmask replace bet with nonlinear registration
          --verbose be verbose
          --qc produce qc images
          --manual-xfm <xfm> use this xfm as initialization
          --search <n> - search radius for patch based segmentation
          --reference <subject id prev> - use this for longitudinal contraint
          --nonbrainconstraint - use non-brain constraint for longitudinal co-registration
         ]
HELP


die $HELP if $#ARGV<1 || !$model_dir;

my $model="$model_dir/model_t1w.mnc";
my $model_bbox="$model_dir/model_t1w_bbox.mnc";
my $model_mask="$model_dir/model_brain.mnc";
my $model_icc ="$model_dir/model_icc_2c.mnc";
my $model_ven ="$model_dir/sample.lst";
my $ven_mask  ="$model_dir/model_avg_vent_90.mnc";
my $lin_ven_mask="$model_dir/model_lin_vent_mask.mnc";
my $nl_ven_mask="$model_dir/model_ven_mask_bbox.mnc";

check_presence($model,$model_mask,$model_ven,$ven_mask,$lin_ven_mask,$nl_ven_mask,$model_bbox);

my ($in,$prefix)=@ARGV;
my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );
do_cmd('mkdir','-p',$prefix,"$prefix/stx","$prefix/clp","$prefix/nu","$prefix/nl","$prefix/vent","$prefix/qc","$prefix/vol");

unless($subject) #figure out subject id
{
  $subject=basename($in,'.gz');
  $subject =~ s/\.gz$//;
  $subject =~ s/\.mnc$//;
}

my $compress=$ENV{MINC_COMPRESS};


#normalize intensity
if($nuc) 
{
  do_cmd('nu_correct',"-iter", 100, "-stop", 0.0001, "-fwhm", 0.1,$in,"$prefix/nu/${subject}.mnc")
   unless -e "$prefix/nu/${subject}.mnc";
  $in="$prefix/nu/${subject}.mnc";
}

if($nuyl)
{
  do_cmd('volume_nuyl',$in,$model,'--fix_zero_padding',"$prefix/clp/${subject}.mnc")
    unless -e "$prefix/clp/${subject}.mnc";
} else {
  do_cmd('volume_pol',$in,$model,'--order',1,'--min',0,'--max',100,'--expfile',"$prefix/clp/${subject}.exp",'--noclamp') 
    unless -e "$prefix/clp/${subject}.exp";

  do_cmd('minccalc','-expfile',"$prefix/clp/${subject}.exp",$in,"$prefix/clp/${subject}.mnc") 
    unless -e "$prefix/clp/${subject}.mnc";
}

unless( -e "$prefix/stx/${subject}.xfm")
{
  #stx registration
  my @args=('bestlinreg.pl',"$prefix/clp/${subject}.mnc",$model,"$tmpdir/stx_${subject}.xfm");
  push @args,'-init_xfm',$manual_xfm if $manual_xfm;
  do_cmd(@args);
  
  # resampling, may be replaced with mincresample --sinc   
  do_cmd('itk_resample','--order',4,"$prefix/clp/${subject}.mnc",
         '--like',$model,'--transform',"$tmpdir/stx_${subject}.xfm",
          "$tmpdir/stx_${subject}.mnc");
  
  #brain extraction , whole head could be used
  if($nlmask)
  {
      do_cmd('nlfit_s','-level',8,"$tmpdir/stx_${subject}.mnc",$model,
        "$tmpdir/nl8_${subject}.xfm");
      
      do_cmd('xfmconcat',"$tmpdir/stx_${subject}.xfm","$tmpdir/nl8_${subject}.xfm","$tmpdir/stx_nl8_${subject}.xfm");
  
      do_cmd('mincresample','-nearest',$model_icc,'-transform',"$tmpdir/stx_nl8_${subject}.xfm",
            "$tmpdir/native_${subject}_mask.mnc",'-invert_transform',
            '-like',"$prefix/clp/${subject}.mnc");
  } else {
    do_cmd('imp_bet.pl',"$tmpdir/stx_${subject}.mnc","$tmpdir/${subject}_mask.mnc");

    do_cmd('mincresample','-nearest',"$tmpdir/${subject}_mask.mnc",
           '-transform',"$tmpdir/stx_${subject}.xfm",
           '-invert_transform',
           "$tmpdir/native_${subject}_mask.mnc",
           '-like',"$prefix/clp/${subject}.mnc");
  }

  # do another linear registration step 
  do_cmd('bestlinreg.pl',
         "$prefix/clp/${subject}.mnc",$model,
         '-init_xfm',"$tmpdir/stx_${subject}.xfm",
         '-source_mask',"$tmpdir/native_${subject}_mask.mnc",
         '-target_mask',$model_icc,
         "$prefix/stx/${subject}.xfm");

  # resampling, may be replaced with mincresample --sinc   
  do_cmd('itk_resample',
         '--order',4,
         "$prefix/clp/${subject}.mnc",
         '--like',$model,
         '--transform',"$prefix/stx/${subject}.xfm",
         "$prefix/stx/${subject}.mnc");

}

if(! -e "$prefix/stx/${subject}_mask.mnc" )
{
  if( $nlmask )
  {
    do_cmd('nlfit_s','-level',8,
            "$prefix/stx/${subject}.mnc",$model,
            "$tmpdir/nl8_2_${subject}.xfm");
            
    do_cmd('mincresample',
          '-nearest',$model_icc,
          '-transform',"$tmpdir/nl8_2_${subject}.xfm",
          "$prefix/stx/${subject}_mask.mnc",
          '-invert_transform',
          '-like',$model_mask);
  } else {
    do_cmd('imp_bet.pl',"$prefix/stx/${subject}.mnc","$prefix/stx/${subject}_mask.mnc");
  }
}


if($whole_head)
{
  # nonlinear fit around ventricles
  do_cmd('nlfit_s','-level',2,
    "$prefix/stx/${subject}.mnc",
    $model,
    '-source_mask',"$prefix/stx/${subject}_mask.mnc",
    '-target_mask',$model_icc,
    "$prefix/nl/${subject}.xfm")
      unless -e "$prefix/nl/${subject}.xfm";
} else {
  # nonlinear fit around ventricles
  do_cmd('nlfit_s','-level',2,
    "$prefix/stx/${subject}.mnc",
    $model,
    '-source_mask',$lin_ven_mask,#"$prefix/stx/${subject}_mask.mnc",
    '-target_mask',$lin_ven_mask,#$model_mask,
    "$prefix/nl/${subject}.xfm")
    unless -e "$prefix/nl/${subject}.xfm";
} 

# search database and build a subject specific average
#do_cmd("$mydir/create_specific_model.rb",
#  "$prefix/stx/${subject}.mnc",
#  "$prefix/stx/${subject}.xfm",
#  "$prefix/nl/${subject}.xfm",
#  $in,
#  "$prefix/vent_/${subject}.mnc",
#  '--model-base',dirname($model),
#  '--model-file',$model_ven,
#  '--model-mask',basename($ven_mask),
#  '--sample',20)
#  unless -e "$prefix/vent_/${subject}.mnc";

unless( -e "$prefix/vent/${subject}.mnc" )
{
  delete $ENV{MINC_COMPRESS} if $compress;

  do_cmd("itk_resample", "$prefix/stx/${subject}.mnc","$tmpdir/nl_${subject}.mnc",
         '--like',$nl_ven_mask,'--transform',"$prefix/nl/${subject}.xfm");

  #renormalize intensities
  do_cmd('volume_pol','--order',1,"$tmpdir/nl_${subject}.mnc",$model_bbox,'--expfile',"$tmpdir/exp",'--noclamp');

  do_cmd('minccalc','-expfile',"$tmpdir/exp","$tmpdir/nl_${subject}.mnc","$tmpdir/nl_${subject}_n.mnc");

  if($reference) # we have got previous timepoint
  {
    do_cmd("itk_resample", "$prefix/stx/${reference}.mnc","$tmpdir/nl_${reference}.mnc",
          '--like',$nl_ven_mask,'--transform',"$prefix/nl/${reference}.xfm");
  
    #renormalize intensities
    do_cmd('volume_pol','--order',1,"$tmpdir/nl_${reference}.mnc",$model_bbox,'--expfile',"$tmpdir/exp_reference",'--noclamp');
    do_cmd('minccalc','-expfile',"$tmpdir/exp_reference","$tmpdir/nl_${reference}.mnc","$tmpdir/nl_${reference}_n.mnc");

    #transform previous ventricle segmentation
    do_cmd('mincresample',"$prefix/stx/${reference}_vent_fuzzy.mnc",
                          '-transform',"$prefix/nl/${reference}.xfm",
                          '-like',$nl_ven_mask,"$tmpdir/nl_vent_${reference}.mnc");
    
    open ADDON,">$tmpdir/previous.txt" or die;
    print ADDON "nl_${reference}_n.mnc,nl_vent_${reference}.mnc";
    close ADDON;

    do_cmd('volume_patches',"$tmpdir/nl_${subject}_n.mnc",
         '--train',$model_ven,"$tmpdir/nl_vent_${subject}.mnc",
         '--mask',$nl_ven_mask,'--search',$search,
         '--patch',1,'--verbose',
         '--train2',"$tmpdir/previous.txt");
  } else {
   
    do_cmd('volume_patches',"$tmpdir/nl_${subject}_n.mnc",
         '--train',$model_ven,"$tmpdir/nl_vent_${subject}.mnc",
         '--mask',$nl_ven_mask,'--search',$search,
         '--patch',1,'--verbose');
  }

  do_cmd('xfmconcat',"$prefix/stx/${subject}.xfm","$prefix/nl/${subject}.xfm","$tmpdir/full_${subject}.xfm");

  $ENV{MINC_COMPRESS}=$compress if $compress;
  
  do_cmd('mincresample',
         '-like',$in,
          "$tmpdir/nl_vent_${subject}.mnc",
         '-invert_transform','-transform',"$tmpdir/full_${subject}.xfm",
         "$tmpdir/vent_${subject}.mnc");
         
  do_cmd('cp',"$tmpdir/vent_${subject}.mnc","$prefix/vent/${subject}_fuzzy.mnc");
  
  do_cmd('mincresample',
         '-like',"$prefix/stx/${subject}.mnc",
         "$tmpdir/nl_vent_${subject}.mnc",
         '-invert_transform','-transform',"$prefix/nl/${subject}.xfm",
         "$tmpdir/stx_vent_${subject}.mnc");
  
  do_cmd('minccalc','-byte','-express','A[0]>0.5?1:0',"$tmpdir/vent_${subject}.mnc","$prefix/vent/${subject}.mnc");

  unless($cleanup || -e "$prefix/stx/${subject}_vent.mnc")
  {
    do_cmd('cp',"$tmpdir/stx_vent_${subject}.mnc","$prefix/stx/${subject}_vent_fuzzy.mnc");
    do_cmd('minccalc','-byte','-express','A[0]>0.5?1:0',"$tmpdir/stx_vent_${subject}.mnc","$prefix/stx/${subject}_vent.mnc");
  }
}

#TODO: finish non-brain contraint scaling factor calculation
if($reference && $nonbrainconstraint&& !(-e "$prefix/clp/lsq12_skull_${subject}_${reference}.xfm" ))
{
  #transform non-brain-masks into native space
  unless( -e "$prefix/clp/lsq12_skull_${subject}_${reference}.xfm")
  {
    extract_skull("$prefix/stx/${subject}.mnc","$prefix/stx/${subject}_skull.mnc") if ! -e "$prefix/stx/${subject}_skull.mnc"; 
    extract_skull("$prefix/stx/${reference}.mnc","$prefix/stx/${reference}_skull.mnc") if ! -e "$prefix/stx/${reference}_skull.mnc"; 
    

    do_cmd('mincresample','-nearest',"$prefix/stx/${subject}_skull.mnc","$tmpdir/native_${subject}_skull.mnc",
          '-transform',"$prefix/stx/${subject}.xfm",'-invert_transform','-like',"$prefix/clp/${subject}.mnc");
    do_cmd('mincresample','-nearest',"$prefix/stx/${reference}_skull.mnc","$tmpdir/native_${reference}_skull.mnc",
          '-transform',"$prefix/stx/${reference}.xfm",'-invert_transform','-like',"$prefix/clp/${reference}.mnc");
  
    #following Jackeline's steps 
    #1  perform lsq12 registration on skull
    do_cmd('bestlinreg_s','-lsq12',"$prefix/clp/${subject}.mnc","$prefix/clp/${reference}.mnc",
          "$tmpdir/lsq12_${subject}_${reference}.xfm");

    #2 perform non-brain part registration, using initialization from head
 
    do_cmd('bestlinreg_s','-lsq12',
          '-init_xfm',"$tmpdir/lsq12_${subject}_${reference}.xfm",
          "$prefix/clp/${subject}.mnc","$prefix/clp/${reference}.mnc",
          '-source_mask',"$prefix/stx/${subject}_skull.mnc",
          '-target_mask',"$prefix/stx/${reference}_skull.mnc",
          "$tmpdir/lsq12_skull_${subject}_${reference}.xfm",'-close');


    # 3 perform backward registration
    do_cmd('bestlinreg_s','-lsq12',"$prefix/clp/${reference}.mnc","$prefix/clp/${subject}.mnc",
          "$tmpdir/lsq12_${reference}_${subject}.xfm");

    #2 perform backward non-brain part registration, using initialization from head
 
    do_cmd('bestlinreg_s','-lsq12',
          '-init_xfm',"$tmpdir/lsq12_${reference}_${subject}.xfm",
          "$prefix/clp/${reference}.mnc","$prefix/clp/${subject}.mnc",
          '-source_mask',"$prefix/stx/${reference}_skull.mnc",
          '-target_mask',"$prefix/stx/${subject}_skull.mnc",
          "$tmpdir/lsq12_skull_${reference}_${subject}.xfm",'-close');

    do_cmd('xfminvert',"$tmpdir/lsq12_skull_${reference}_${subject}.xfm","$tmpdir/lsq12_skull_${reference}_${subject}_inv.xfm");
    do_cmd('xfmavg',"$tmpdir/lsq12_skull_${subject}_${reference}.xfm","$tmpdir/lsq12_skull_${reference}_${subject}_inv.xfm",
           "$prefix/clp/lsq12_skull_${subject}_${reference}.xfm");

  }

}

if($reference && $nonbrainconstraint)
{
  my @ts=split(/\s+/,`xfm2param $prefix/clp/lsq12_skull_${subject}_${reference}.xfm|fgrep scale`);
  $corr=$ts[1]*$ts[2]*$ts[3];
  
  my @tr=split(/\s+/,`xfm2param $prefix/clp/lsq12_skull_${subject}_${reference}.xfm|fgrep translation`);
  $translation=sqrt($tr[1]*$tr[1]+$tr[2]*$tr[2]+$tr[3]*$tr[3]);
}

my @sc=split(/\s+/,`xfm2param $prefix/stx/${subject}.xfm|fgrep scale`);
my $scale=$sc[1]*$sc[2]*$sc[3];

#calculate ventricle volume
my $lvv=`mincstats -q -sum $prefix/vent/${subject}_fuzzy.mnc`;
$lvv*=1.0;
my @vox=split("\n",`mincinfo -attvalue xspace:step -attvalue yspace:step -attvalue zspace:step  $prefix/vent/${subject}_fuzzy.mnc`);
$lvv*=$vox[0]*$vox[1]*$vox[2];

open OUT,">$prefix/vol/${subject}_lvv.csv" or die;
print OUT "${subject},${lvv},$scale,$corr,$translation\n";
close OUT;


if($qc)
{
  if($whole_head)
  {
  #"$prefix/stx/${subject}_mask.mnc"
  
    do_cmd('minccalc', '-byte', '-express','clamp(A[0]/200,0,1)',
          "$prefix/stx/${subject}.mnc","$tmpdir/${subject}_stx.mnc");
    
    do_cmd('minccalc', '-express','A[1]>0.5?0:A[0]',"$prefix/stx/${subject}_mask.mnc",
          "$prefix/stx/${subject}_vent.mnc", "$tmpdir/${subject}_mod.mnc");
  
    do_cmd('minc_qc_rgb.pl',"$tmpdir/${subject}_mod.mnc", "$tmpdir/${subject}_stx.mnc", 
          "$prefix/stx/${subject}_vent.mnc","$prefix/qc/${subject}.jpg" ) 
          unless -e "$prefix/qc/${subject}.jpg";
    
  } else {
    do_cmd('minc_qc.pl',"$prefix/stx/${subject}.mnc",'--image-range',0,100,
          '--mask',"$prefix/stx/${subject}_vent.mnc","$prefix/qc/${subject}.jpg" ) 
          unless -e "$prefix/qc/${subject}.jpg";
  }
  
  if($reference && $nonbrainconstraint)
  {
    do_cmd('mincresample',"$prefix/clp/${subject}.mnc","$tmpdir/${subject}_${reference}_align.mnc",
           '-like',"$prefix/clp/${reference}.mnc",
           '-transform',"$prefix/clp/lsq12_skull_${subject}_${reference}.xfm");
    do_cmd('minc_qc.pl','--image-range',0,100,'--red',
           '--green-mask','--mask-range',0,100,"$prefix/clp/${reference}.mnc",
           '--mask',"$tmpdir/${subject}_${reference}_align.mnc",'--big',
           "$prefix/qc/ref_${subject}_${reference}.jpg") unless -e "$prefix/qc/ref_${subject}_${reference}.jpg";
  }
}

if($cleanup)
{
  unlink( "$prefix/nl/${subject}_4mm.xfm",
          "$prefix/nl/${subject}.xfm",
          "$prefix/nl/${subject}_4mm_grid_0.mnc",
          "$prefix/nl/${subject}_grid_0.mnc",
          "$prefix/stx/${subject}_mask.mnc",
          "$prefix/clp/${subject}.mnc",
          "$prefix/nu/${subject}.mnc",
          "$prefix/stx/${subject}_vent.mnc" );
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

sub extract_skull
{
  my ($in,$out)=@_;
  do_cmd('mincbet',$in,"$tmpdir/bet",'-s','-m','-n');
  do_cmd('itk_morph','--exp','D[5] D[5]', "$tmpdir/bet_mask.mnc","$tmpdir/bet_mask_d.mnc",'--clob');
  do_cmd('itk_morph','--exp', 'D[2]', "$tmpdir/bet_skull.mnc","$tmpdir/bet_skull_d.mnc",'--clob');
  do_cmd('minccalc','-byte','-express','A[0]>0.5&&A[1]>0.5&&A[2]<0.5?1:0',
         "$tmpdir/bet_skull_d.mnc","$tmpdir/bet_mask_d.mnc","$tmpdir/bet_mask.mnc","$tmpdir/all_skull.mnc",'-clob');
  do_cmd('mincreshape','-dimrange','zspace=60,120',"$tmpdir/all_skull.mnc","$tmpdir/skull.mnc",'-clob');
  do_cmd('mincresample','-nearest', "$tmpdir/skull.mnc",'-like',$in,'-byte','-clob',$out);
  do_cmd('rm','-f',"$tmpdir/skull.mnc","$tmpdir/all_skull.mnc","$tmpdir/bet_skull_d.mnc","$tmpdir/bet_mask_d.mnc","$tmpdir/bet_mask.mnc");
}

# kate: space-indent on; indent-width 2; replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl perl
