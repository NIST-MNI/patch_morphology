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
my $qc;
my $manual_xfm;
my $search=2;
my $patch=1;
my $exclude=0;
my $classes=18;
my $threshold=0.0;
my $variant='labels';
my $compare;

GetOptions (
          'verbose'        => \$verbose,
          'qc'             => \$qc,
          'clobber'        => \$clobber,
          'model-dir=s'    => \$model_dir,
          'subject=s'      => \$subject,
          'search=n'       => \$search,
          'patch=n'        => \$patch,
          'exclude'        => \$exclude,
          'classes=n'      => \$classes,
          'variant=s'      => \$variant,
          'threshold=f'    => \$threshold,
          'compare=s'      => \$compare,
          );

my $HELP=<<HELP ;
Reference: Pierrick Coupé, , José V. Manjón, Vladimir Fonov, Jens Pruessner, Montserrat Robles, D. Louis Collins 
           "Patch-based segmentation using expert priors: Application to hippocampus and ventricle segmentation" 
           http://dx.doi.org/10.1016/j.neuroimage.2010.09.018

Usage $me input_t1.mnc output_prefix  
          --model-dir <anatomical model directory> - required!
         [
          --search <n> - search radius for patch based segmentation
          --patch <n> - patch radius
          --subject <id> - subject id
          --qc - produce qc images
          --exclude exclude subject from the library for cross-validation, WARNING: will create a temporary file in library directory
          --variant <s>
          --compare <s>
         ]
HELP


die $HELP if $#ARGV<1 || !$model_dir;

my $model     ="$model_dir/model_t1w.mnc";
my $model_mask="$model_dir/model_t1w_mask.mnc";
my $model_bbox="$model_dir/model_t1w_bbox.mnc";
my $model_bbox_mask="$model_dir/model_t1w_bbox_mask.mnc";
my $model_outline="$model_dir/model_t1w_outline.mnc";


my $model_samples = "$model_dir/samples.lst";
my $model_samples_ex = "$model_dir/samples_${subject}.lst";

check_presence($model,$model_mask,$model_samples,$model_bbox,$model_bbox_mask);

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


#linear registration
if(! -e "${prefix}/stx_${subject}.xfm")
{
  do_cmd('bestlinreg_s',
         $in,$model,'-nmi',
         "${prefix}/stx_${subject}.xfm");
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
  do_cmd('minc_qc.pl',"${prefix}/stx_${subject}.mnc",'--mask',$model_outline,"${prefix}/qc_stx_${subject}.jpg");
}

# preliminary non-linear registration to get approximate mask
if(! -e "${prefix}/stx_${subject}_mask.mnc")
{
  do_cmd('nlfit_s',
         '-level',4,'-start',16,
         "${prefix}/stx_${subject}.mnc",$model,
         "$tmpdir/rough.xfm");
  
  do_cmd('itk_resample','--labels',$model_bbox_mask,'--transform',"$tmpdir/rough.xfm",'--invert_transform',
         "${prefix}/stx_${subject}_mask.mnc",'--byte','--like',$model);
  
}


# second intensity normalization

if(! -e "${prefix}/stx_${subject}_bbox.mnc")
{
  do_cmd('volume_pol',"${prefix}/stx_${subject}.mnc",$model_bbox,
         '--source_mask',"${prefix}/stx_${subject}_mask.mnc",
         '--target_mask',$model_bbox_mask,
         '--order',0,'--expfile',"$tmpdir/normalize.txt",'--noclamp');
  
  do_cmd('minccalc','-expfile',"$tmpdir/normalize.txt","${prefix}/stx_${subject}.mnc","${tmpdir}/stx_${subject}_n.mnc");
  do_cmd('mincresample',"${tmpdir}/stx_${subject}_n.mnc","${prefix}/stx_${subject}_bbox.mnc",'-nearest','-like',$model_bbox);
}


# non-linear registration
if(! -e "${prefix}/nl_${subject}.xfm")
{
  
  do_cmd('nlfit_s',
         '-level',2,'-start',16,"${prefix}/stx_${subject}_bbox.mnc",$model_bbox,
         "${prefix}/nl_${subject}.xfm");
}

if(! -e "${prefix}/nl_full_${subject}.xfm")
{
  do_cmd('xfmconcat',"${prefix}/stx_${subject}.xfm","${prefix}/nl_${subject}.xfm","${prefix}/nl_full_${subject}.xfm");
}

if(! -e "${prefix}/nl_${subject}.mnc")
{
  do_cmd('itk_resample',
         $in,'--like',$model_bbox_mask,
         '--transform',"${prefix}/nl_full_${subject}.xfm",
         "${prefix}/nl_${subject}.mnc" );
}

if(! -e "${prefix}/nl_${subject}_n.mnc")
{
  do_cmd('volume_pol',"${prefix}/nl_${subject}.mnc",$model_bbox,
         '--source_mask',$model_bbox_mask,'--target_mask',$model_bbox_mask,
         '--order',0,'--expfile',"$tmpdir/normalize2.txt",'--noclamp');
  
  do_cmd('minccalc','-expfile',"$tmpdir/normalize2.txt","${prefix}/nl_${subject}.mnc","${prefix}/nl_${subject}_n.mnc");
}

if(! -e "$prefix/nl_${subject}_${variant}.mnc" )
{
  my $train=$model_samples;
  if($exclude)
  {
    do_cmd("fgrep -v $subject $model_samples>$model_samples_ex");
    $train=$model_samples_ex;
  }
  my @args=('itk_patch_morphology',
          "${prefix}/nl_${subject}_n.mnc",
          '--mask',$model_bbox_mask,
          '--train',$train,
          '--search',$search,
          '--patch',$patch,'--discrete',$classes,
          "$prefix/nl_${subject}_${variant}.mnc");
  push(@args,'--threshold',$threshold)  if $threshold>0.0;
  
  do_cmd(@args);

  if($exclude)
  {
    do_cmd('rm','-f',$model_samples_ex);
  }
}

# transform labels back into original space...
if(! -e "$prefix/${subject}_${variant}.mnc" )
{
  do_cmd('itk_resample','--labels',"$prefix/nl_${subject}_${variant}.mnc","$prefix/${subject}_${variant}.mnc",'--like',$in,
         '--transform',"${prefix}/nl_full_${subject}.xfm",'--invert_transform','--byte');
  
  #do_cmd('mincreshape','-byte','-image_range',0,255,'-valid_range',0,255,"$tmpdir/${subject}_${variant}.mnc","$prefix/${subject}_${variant}.mnc");
}

if($qc && ! -e "${prefix}/qc_${subject}_labels.jpg")
{
  do_cmd('minc_qc.pl',$in,'--mask',"$prefix/${subject}_${variant}.mnc","${prefix}/qc_${subject}_labels.jpg");
}

if($compare && (! -e "$prefix/compare_${subject}_${variant}.txt"))
{
  my $m=`volume_gtc_similarity --csv $prefix/${subject}_${variant}.mnc $compare `;
  chomp($m);
  open OUT,">$prefix/compare_${subject}_${variant}.txt" or die "Can't open $compare for writing\n";
  print OUT "$subject,$search,$patch,$threshold,$m\n";
  close OUT;
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

# kate: space-indent on; indent-width 2; replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl perl
