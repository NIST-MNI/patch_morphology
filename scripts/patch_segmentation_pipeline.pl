#! /usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Parallel::ForkManager;

my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $mydir=dirname($0);
my $model_dir;
my $nuc;
my $cleanup=0;
my $subject;
my $qc;
my $manual_xfm;

my $search=2;
my $patch=1;
my $exclude=0;
my $classes=36;
my $threshold=0.0;
my $variant='labels';
my $compare;
my $short;
my $preprocess=0;
my $build;
my $symmetric;
my $preselect=0;

my $use_minctracc=0;
my $use_ants=0;
my $resample_labels_order=1;
my $pairwise;
my $refine;
my $qc_lut;
my $nuc;
my $nuc3T;
my $baa=0;
my $no_nuyl=0;
my $no_gco=0;
my $no_nl=0;
my $list='samples.lst';
my $nomask=0;
my $patch_normalize=0;

GetOptions (
          'verbose'        => \$verbose,
          'qc=s'           => \$qc,
          'clobber'        => \$clobber,
          'cleanup'        => \$cleanup,
          'model-dir=s'    => \$model_dir,
          'subject=s'      => \$subject,
          'search=n'       => \$search,
          'patch=n'        => \$patch,
          'exclude'        => \$exclude,
          'classes=n'      => \$classes,
          'variant=s'      => \$variant,
          'threshold=f'    => \$threshold,
          'compare=s'      => \$compare,
          'short'          => \$short,
          'preprocess'     => \$preprocess,
          'minctracc'      => \$use_minctracc,
          'ants'           => \$use_ants,
          'build=s'        => \$build,
          'symmetric=s'    => \$symmetric,
          'resample_labels_order=n' =>\$resample_labels_order,
          'preselect=n'    => \$preselect,
          'pairwise'       => \$pairwise,
          'refine'         => \$refine,
          'qc_lut=s'       => \$qc_lut,
          'nuc'            => \$nuc,
          'nuc3T'          => \$nuc3T,
          'no_nuyl'        => \$no_nuyl,
          'no_gco'         => \$no_gco,
          'no_nl'          => \$no_nl,
          'baa'            => \$baa,
          'list=s'         => \$list,
          'nomask'         => \$nomask,
          'patch_normalize'=> \$patch_normalize, 
          );

my $HELP=<<HELP ;
Reference: Katrin Weier, Vladimir Fonov, Karyne Lavoie, Julien Doyon and D. Louis Collins
          "Rapid automatic segmentation of the human cerebellum and its lobules (RASCAL) 
          Implementation and application of the patch-based label-fusion technique with a template library to segment the human cerebellum"
          http://dx.doi.org/10.1002/hbm.22529
          
Usage $me input_t1.mnc output_prefix  
          --model-dir <anatomical model directory> - required!
         [
          --search <n> - search radius for patch based segmentation
          --patch <n> - patch radius
          --subject <id> - subject id
          --preselect <n> - number of samples to preselect, default: 0 - all
          --qc_lut <file.lut> lut file for qc image
          --exclude exclude subject from the library for cross-validation, WARNING: will create a temporary file in library directory
          --variant <s>
          --compare <s>
          --short - assume labels are 16 bit
          --preprocess - perform stx registration first
          --nuc - run N3 in preprocessing stage, for 1.5T
          --nuc3T - run N3 in preprocessing stage, for 3T
          --resample_labels_order <n> use Nth order for resampling labels
          --baa use anti-aliasing while resampling labels
          --no_nuyl don't run NUYL normalization
          --no_gco  disable graph-cut smoothing
          --cleanup  remove unneeded temporary files
          --list <list name> use alternative list name instead of samples.lst
          --nomask - disable masking for non-linear registration and patch application
          --patch_normalize apply patch-based normalization instead of the usual
          =====================
          Nonlinear registration options: by default perform only linear registration !
          --ants - use mincANTS for nonlinear normalization
          --mintracc - use minctracc for nonlinear normalization
          --elastix - use elastix [TODO]
          --pairwise - perform pairwize coregistration 
          --refine - instead of applying existing xfm refine it!
          =====================
          Library building options
          --build <sample>  - instead of performing segmentation, use recovered transformations to build the library
          --symmetric <remap.lut> - assume symmetric model, remap both sample and t1
         ]
HELP


die $HELP if $#ARGV<1 || !$model_dir;

my $model          ="$model_dir/model_t1w.mnc";
my $model_mask     ="$model_dir/model_mask.mnc";

my $model_bbox="$model_dir/model_t1w_bbox.mnc";
my $model_bbox_mask="$model_dir/model_mask_bbox.mnc";

my $model_samples = "$model_dir/$list";
my $model_samples_ex = "$model_dir/samples_${subject}.lst";

my $model_remap_fwd= "$model_dir/remap_fwd.lut";
my $model_remap_back= "$model_dir/remap_back.lut";
my $model_label_energy="$model_dir/label_interaction.csv";
my $model_patch_db="$model_dir/samples_reduced.db";
my $model_patch_index="$model_dir/samples_reduced.index";
my $remap=0;

check_presence($model,$model_mask,$model_bbox,$model_bbox_mask,$model_label_energy);

check_presence($model_samples) unless $build; # we don't need existing samples if we are building the library

$remap=1 if  -e $model_remap_fwd && -e $model_remap_back; # assume that we are going to be mapping labels forward and backward!


my ($in,$prefix)=@ARGV;
my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );
do_cmd('mkdir','-p',$prefix);
do_cmd('mkdir','-p',"$prefix/${subject}/");

unless($subject) #figure out subject id
{
  $subject=basename($in,'.gz');
  $subject =~ s/\.gz$//;
  $subject =~ s/\.mnc$//;
}

my $compress=$ENV{MINC_COMPRESS};
#$ENV{ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}=1;
my $in_original=$in;

#linear registration
if($preprocess)
{
  if($nuc)
  {
    if(! -e "${prefix}/nuc_${subject}.mnc")
    {
      do_cmd('pipeline_correct.pl', $in,
             "${prefix}/nuc_${subject}.mnc",
             '--model', $model,
             '--model-mask', $model_mask);
    }
    $in="${prefix}/nuc_${subject}.mnc";
  } elsif($nuc3T) {
    if(! -e "${prefix}/nuc_${subject}.mnc")
    {
      do_cmd('pipeline_correct.pl', $in,
             "${prefix}/nuc_${subject}.mnc",
             '--model', $model,
             '--model-mask', $model_mask,
             '--3t');
    }
    $in="${prefix}/nuc_${subject}.mnc";
  }
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

  # preliminary non-linear registration to get approximate mask
  if(! -e "${prefix}/stx_${subject}_mask.mnc")
  {
    do_cmd('nlfit_s',
          '-level',4,'-start',16,
          "${prefix}/stx_${subject}.mnc",$model,
          "$tmpdir/rough.xfm");
    
    do_cmd('itk_resample','--labels',$model_mask,'--transform',"$tmpdir/rough.xfm",'--invert_transform',
          "${prefix}/stx_${subject}_mask.mnc",'--byte','--like',$model);
    
  }

  # second intensity normalization
  if(! -e "${prefix}/stx_${subject}_n.mnc")
  {
      do_cmd('volume_pol',"${prefix}/stx_${subject}.mnc",$model,
              '--source_mask',"${prefix}/stx_${subject}_mask.mnc",
              '--target_mask',$model_mask,
              '--order',0,'--expfile',"$tmpdir/normalize.txt",
              '--noclamp');

      do_cmd('minccalc','-expfile',"$tmpdir/normalize.txt","${prefix}/stx_${subject}.mnc","${prefix}/stx_${subject}_n.mnc");
  }
  $in="${prefix}/stx_${subject}_n.mnc";
}

if(! -e "${prefix}/stx_bbox_${subject}.mnc")
{
  do_cmd('itk_resample',
          $in,'--like',$model_bbox_mask,
          "${prefix}/stx_bbox_${subject}.mnc");
}

if(! -e "${prefix}/stx_bbox_${subject}_n.mnc")
{
  do_cmd('itk_minc_nonlocal_filter','--beta',0.5,'--anlm',
          "${prefix}/stx_bbox_${subject}.mnc",
          "${tmpdir}/stx_bbox_${subject}_den.mnc");

  if ($no_nuyl)
  {
   do_cmd('cp',
            "${tmpdir}/stx_bbox_${subject}_den.mnc",
            "${prefix}/stx_bbox_${subject}_n.mnc");
  } else {
    if ($patch_normalize)
    {
      do_cmd('flann_patch_normalize.pl', 
              '--mask', $model_bbox_mask, 
              '--db', $model_patch_db, 
              '--index', $model_patch_index, 
              "${tmpdir}/stx_bbox_${subject}_den.mnc", 
              "${prefix}/stx_bbox_${subject}_n.mnc");
    } else {
      do_cmd('minc_nuyl',
              "${tmpdir}/stx_bbox_${subject}_den.mnc",
              $model_bbox,
              "${prefix}/stx_bbox_${subject}_n.mnc",
              '--source-mask',$model_bbox_mask,
              '--target-mask',$model_bbox_mask,
              '--linear'
            );
    }
  }
}

if($build && $symmetric)
{
  if(!-e "${tmpdir}/flip_${subject}_${variant}.xfm")
  {
    do_cmd('param2xfm','-scale',-1,1,1,"${tmpdir}/flip_${subject}_${variant}.xfm");
  }

  if(!-e "${prefix}/stx_bbox_${subject}_n_flip.mnc")
  {
    do_cmd('itk_resample',
          "${prefix}/stx_bbox_${subject}_n.mnc",
          "${prefix}/stx_bbox_${subject}_n_flip.mnc",
          '--transform',"${tmpdir}/flip_${subject}_${variant}.xfm",'--order',0);
  }
  
  unless( -e "${prefix}/nl_${subject}_flip.xfm" || $no_nl ) {
    nl_reg("${prefix}/stx_bbox_${subject}_n_flip.mnc",$model_bbox,$model_bbox_mask,"${prefix}/nl_${subject}_flip.xfm");
  }
}

unless($no_nl)
{
 unless(-e "${prefix}/nl_${subject}.xfm" || $pairwise ) {
   nl_reg("${prefix}/stx_bbox_${subject}_n.mnc",$model_bbox,$model_bbox_mask,"${prefix}/nl_${subject}.xfm");
   }
}

if($build)
{
  my @args=('itk_resample',
            '--labels',$build,
            "${prefix}/${subject}_labels.mnc",
            '--like',$model_bbox_mask,
            );

  push (@args,'--lut',$model_remap_fwd) if $remap;

  if($short) {
    push (@args,'--short');
  } else {
    push (@args,'--byte');
  }
  do_cmd(@args) unless -e "${prefix}/${subject}_labels.mnc";

  if($symmetric) # flip the sample and the labels
  {
    if(!-e "${prefix}/${subject}_labels_flip.mnc")
    {
      @args=('itk_resample',
           "${prefix}/${subject}_labels.mnc",
           "${prefix}/${subject}_labels_flip.mnc",
           '--transform',"$tmpdir/flip_${subject}_${variant}.xfm",
           '--lut',$symmetric,'--labels');
      
      if($short) {
        push (@args,'--short');
      } else {
        push (@args,'--byte');
      }
      do_cmd(@args);
    }
  }
} else {
  if(! -e "${prefix}/stx_bbox_${subject}_${variant}_dumb.mnc" )
  {
    my $train=$model_samples;
    
    if($exclude) #TODO: move down 
    {
      do_cmd("fgrep -v $subject $model_samples>$model_samples_ex");
      $train=$model_samples_ex;
    }
    my $custom_train="${prefix}/${subject}/samples.lst";
    
    open TRAIN,"<$train" or die "Can't open $train\n";
    open CTRAIN,">$custom_train" or die "Can't open $custom_train\n";
    my $line;

    my $pm = Parallel::ForkManager->new($ENV{ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS});
    
    foreach $line (<TRAIN>)  
    {   
      chomp($line);
      my ($bbox,$lab,$xfm)=split(/,/,$line);
      my $n=basename($bbox,'.mnc');
      
      my $full_xfm="${tmpdir}/full_${n}_nl.xfm";
      print CTRAIN "${n}.mnc,${n}_lab.mnc\n";

      $pm->start() and next; # do the fork

      if(! -e "${prefix}/${subject}/${n}.mnc" )
      {
        if($pairwise)
        {
          # TODO: use previous xfm to initialize?
          $full_xfm="${prefix}/${subject}/full_${n}_nl.xfm";
          unless(-e $full_xfm) {
            nl_reg("${prefix}/stx_bbox_${subject}_n.mnc","$model_dir/$bbox",$model_bbox_mask,$full_xfm);
          }
        } else {
          # TODO: run xfm_normalize ?
          do_cmd('xfminvert',"$model_dir/$xfm","${tmpdir}/temp1_${n}_nl.xfm");
          do_cmd('xfmconcat',"${prefix}/nl_${subject}.xfm","${tmpdir}/temp1_${n}_nl.xfm","${tmpdir}/temp_${n}_nl.xfm");
          if($refine)
          {
            do_cmd('xfm_normalize.pl','--like',"${prefix}/stx_bbox_${subject}_n.mnc",'--exact',"${tmpdir}/temp_${n}_nl.xfm","${tmpdir}/temp2_${n}_nl.xfm");
            nl_reg_refine("${prefix}/stx_bbox_${subject}_n.mnc","$model_dir/$bbox",$model_bbox_mask,"${tmpdir}/temp2_${n}_nl.xfm",$full_xfm);
            do_cmd('rm','-f',"${tmpdir}/temp_${n}_nl_grid_?.mnc","${tmpdir}/temp2_${n}_nl_grid_?.mnc");
          } else {
            do_cmd('xfm_normalize.pl','--like',"${prefix}/stx_bbox_${subject}_n.mnc",'--exact',"${tmpdir}/temp_${n}_nl.xfm",$full_xfm);
            do_cmd('rm','-f',"${tmpdir}/temp_${n}_nl_grid_?.mnc","${tmpdir}/temp2_${n}_nl_grid_?.mnc");
          }
          do_cmd('rm','-f',"${tmpdir}/temp1_${n}_nl.xfm","${tmpdir}/temp1_${n}_nl_grid_?.mnc");
        }
      
      
        if(! -e "${prefix}/${subject}/${n}.mnc")
        {
          do_cmd('itk_resample',"$model_dir/$bbox",'--like',"${prefix}/stx_bbox_${subject}_n.mnc",
               '--transform',$full_xfm,'--invert_transform',"${prefix}/${subject}/${n}.mnc");
        }

        my @args=('itk_resample',"$model_dir/$lab",'--like',"${prefix}/stx_bbox_${subject}_n.mnc",
                  '--transform',$full_xfm,"${prefix}/${subject}/${n}_lab.mnc",
                  '--labels','--invert_transform');
        
        push (@args,'--baa','--order',2) if $baa;

        if($short){
          push (@args,'--short');
        } else {
          push (@args,'--byte');
        }
        #TODO: order?
      
        do_cmd(@args) if ! -e "${prefix}/${subject}/${n}_lab.mnc";
      }
      
      
      unlink('rm','-f',"${tmpdir}/${n}_nl_grid_?.mnc");# TODO: check this
      $pm->finish();
    }
    $pm->wait_all_children();
    close TRAIN;
    close CTRAIN;
    
    # warp the mask to match the subject
    do_cmd('itk_resample','--byte','--transform',"${prefix}/nl_${subject}.xfm",
                          '--labels','--invert_transform',
                          $model_bbox_mask,"${tmpdir}/mask_${subject}.mnc");
           
    # now create individual library by applying transformations to each sample
    my @args=('itk_patch_morphology',
            "${prefix}/stx_bbox_${subject}_n.mnc",
            '--train',$custom_train,
            '--search',$search,
            '--patch',$patch,'--discrete',$classes,
            '--prob',"${prefix}/${subject}/${subject}_prob", 
            '--mask',"${tmpdir}/mask_${subject}.mnc",
            "$prefix/stx_bbox_${subject}_${variant}_dumb.mnc");

#    push(@args,'--mask',$model_bbox_mask) unless $nomask;
    push(@args,'--threshold',$threshold,'--iter',2)  if $threshold>0.0;
    push(@args,'--top',$preselect) if $preselect>0;

    do_cmd(@args) unless -e "$prefix/stx_bbox_${subject}_${variant}_dumb.mnc";

    if($exclude)
    {
      do_cmd('rm','-f',$model_samples_ex);
    }
    my @probs;
    my $i;
    for ($i=0;$i<$classes;$i+=1)
    {
      push @probs, sprintf("${prefix}/${subject}/${subject}_prob_%02d.mnc", $i);
    }
    
    do_cmd('gco_classify', '--label-energy', $model_label_energy, @probs, 
      "$prefix/stx_bbox_${subject}_${variant}_gco.mnc", 
      '--iter',1000,'--mask',"${tmpdir}/mask_${subject}.mnc");
    
    # cleanup
    do_cmd('rm','-rf',"${prefix}/${subject}") if $cleanup;
  }
  my $seg_result=$no_gco?"$prefix/stx_bbox_${subject}_${variant}_dumb.mnc":"$prefix/stx_bbox_${subject}_${variant}_gco.mnc";
  
  if(! -e "${prefix}/stx_bbox_${subject}_${variant}.mnc" && $remap)
  {
    my @args=('itk_resample',
              '--labels',
              $seg_result,
              "$prefix/stx_bbox_${subject}_${variant}.mnc",
              '--lut', $model_remap_back);

    if($short){
      push (@args,'--short');
    } else {
      push (@args,'--byte');
    }
    do_cmd(@args);
  }

  my $ref=$in;
  $ref=$in_original if($preprocess);
  
  # transform labels back into original space...
  if(! -e "$prefix/${subject}_${variant}.mnc" )
  {
    my @args=('itk_resample',
              '--labels',
              $seg_result,
              "$prefix/${subject}_${variant}.mnc",
              '--like',$ref);
    push (@args,'--baa','--order',2) if $baa;
    push (@args,'--lut',$model_remap_back) if $remap;
    push (@args,'--transform', "${prefix}/stx_${subject}.xfm" ,'--invert') if $preprocess;

    if($short){
      push (@args,'--short');
    } else {
      push (@args,'--byte');
    }
    do_cmd(@args);
  }

  if($compare )
  {
    if(! -e "$prefix/compare_${subject}_${variant}.txt" )
    {
      print STDOUT "volume_gtc_similarity --csv $prefix/${subject}_${variant}.mnc $compare";
      
      my $m=`volume_gtc_similarity --csv $prefix/${subject}_${variant}.mnc $compare`;
      chomp($m);
      open OUT,">$prefix/compare_${subject}_${variant}.txt" or die "Can't open $compare for writing\n";
      print OUT "$subject,$search,$patch,$threshold,$m\n";
      close OUT;
    }
    
    if(! -e "$prefix/stx_bbox_${subject}_${variant}_diff.mnc" )
    {
      my @args=('itk_resample',
                '--labels',$compare,
                "${tmpdir}/compare_${subject}_labels.mnc",
                '--like',$model_bbox_mask,
                );
      push (@args,'--lut',$model_remap_fwd) if $remap;

      if($short) {
        push (@args,'--short');
      } else {
        push (@args,'--byte');
      }
      do_cmd(@args);
      
      do_cmd('minccalc','-byte','-express','abs(A[0]-A[1])<0.5?0:1',
             "${tmpdir}/compare_${subject}_labels.mnc",
             $seg_result, 
             "$prefix/stx_bbox_${subject}_${variant}_diff.mnc");
    }
    
    if(! -e "$prefix/nl_${subject}_${variant}_diff.mnc" )
    {
      do_cmd('itk_resample','--byte',
             "$prefix/stx_bbox_${subject}_${variant}_diff.mnc",
             "$prefix/nl_${subject}_${variant}_diff.mnc",
             '--transform',"${prefix}/nl_${subject}.xfm",'--order',1);
    }
    
    
    do_cmd('minc_qc.pl', "${prefix}/stx_bbox_${subject}_n.mnc",
           '--image-range', 0, 150,
           '--mask',"$prefix/stx_bbox_${subject}_${variant}_diff.mnc",
           '--big','--title',"${subject}_${variant}",
           "$prefix/qc_diff_${subject}_${variant}.jpg")
      unless -e "$prefix/qc_diff_${subject}_${variant}.jpg";
  }

  unless( -e "$prefix/qc_seg_${subject}_${variant}.jpg")
  {
    my $remapped=$seg_result;

    if($remap)
    {
      do_cmd('itk_resample','--labels','--short',
            '--lut',$model_remap_back,
             $seg_result,
             "$tmpdir/stx_bbox_${subject}_${variant}_remap.mnc");
      
      $remapped="$tmpdir/stx_bbox_${subject}_${variant}_remap.mnc";
    }

    my @args=('minc_qc.pl', "${prefix}/stx_bbox_${subject}_n.mnc",
              '--mask', $remapped,
              '--image-range', 0, 150,
              "$prefix/qc_seg_${subject}_${variant}.jpg",
              '--big', '--discrete-mask',
              '--title',"${subject}_${variant}");

    push @args,('--mask-lut', $qc_lut) if $qc_lut ;

    do_cmd(@args) 
  }
}


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    print STDOUT join(',',@_)."\n" if $verbose;
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


sub nl_reg {
  my ($in,$target,$mask,$out)=@_;
  my @args;
  if($use_minctracc)
  {
    # non-linear registration
    @args=('nlfit_s',
            '-level',2,'-start',16,
             $in,
             $target,
             $out
            );
    push(@args,'-source_mask',$mask,'-target_mask',$mask) unless $nomask;

  } elsif($use_ants) {
    #options for high resolution (0.5x0.5x0.5mm voxels)

    @args=('ANTS',3,
            '-i','100x100x50',
            '-m',"CC[$in,$target,1,4]",
            '-r','Gauss[1.0,1.0]',
            '-t','SyN[0.25]',
            '-o',$out);
    push(@args,'-x',$mask)  unless $nomask;
    
  } else {
    print "\n\nSkipping non-linear registration!\n\n";
    @args=('bestlinreg_s','-nmi','-close',
          $in,$target,$out );
    push(@args,'-source_mask',$mask,'-target_mask',$mask) unless $nomask;
  }

  do_cmd(@args);
}

sub nl_reg_refine {
  my ($in,$target,$mask,$xfm,$out)=@_;
  
  do_cmd('nlfit_l',
        '-level',2,'-init_xfm',$xfm,
         $in,
         $target,
         $out,
         '-source_mask',$mask,
         '-target_mask',$mask
        );
}

# kate: space-indent on; indent-width 2; replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2; hl perl
