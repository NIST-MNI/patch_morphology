#! /bin/sh

bindir=$1
rundir=$2
threshold=2

if [ -z $rundir ];then
  echo "Usage $0 bindir rundir"
  exit 1 
fi

BC=$(which bc)

if [ -z $BC ];then
  echo bc is missing! >&2
  exit 1
fi

echo Running patch-based segmentation test
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1

make_phantom -ellipse \
             -continuous -short \
             -start -50 -50 -50 \
             -nelements 50 50 50 \
             -background 0 \
             -step 2 2 2 \
             -fill_value 100 \
             -edge_value 50  \
             -width 10 10 10 \
             -clob \
             -center -30 0 0  \
             $rundir/test_perfect_sphere_1.mnc

param2xfm -shift 30 0 0 $rundir/shift_x30.xfm -clob

mincresample -nearest -transform $rundir/shift_x30.xfm -use_input $rundir/test_perfect_sphere_1.mnc $rundir/test_perfect_sphere_2.mnc -clob
mincresample -nearest -transform $rundir/shift_x30.xfm -use_input $rundir/test_perfect_sphere_2.mnc $rundir/test_perfect_sphere_3.mnc -clob

minccalc -express 'A[0]>1?1:0' -byte -label $rundir/test_perfect_sphere_1.mnc $rundir/test_perfect_sphere_1_A.mnc -clob
minccalc -express 'A[0]>1?2:0' -byte -label $rundir/test_perfect_sphere_2.mnc $rundir/test_perfect_sphere_2_B.mnc -clob
minccalc -express 'A[0]>1?3:0' -byte -label $rundir/test_perfect_sphere_3.mnc $rundir/test_perfect_sphere_3_C.mnc -clob


random_volume --gauss 10 $rundir/test_perfect_sphere_1.mnc $rundir/test_noise_seg.mnc  --clob

# generate segmentation sample with some noise
minccalc -express 'A[0]+A[1]+A[2]+A[3]' $rundir/test_perfect_sphere_1.mnc $rundir/test_perfect_sphere_2.mnc $rundir/test_perfect_sphere_3.mnc $rundir/test_noise_seg.mnc $rundir/test_seg_sample.mnc

#generate training "library"
echo test_perfect_sphere_1.mnc,test_perfect_sphere_1_A.mnc >  $rundir/patch_seg_library.txt
echo test_perfect_sphere_2.mnc,test_perfect_sphere_2_B.mnc >> $rundir/patch_seg_library.txt
echo test_perfect_sphere_3.mnc,test_perfect_sphere_3_C.mnc >> $rundir/patch_seg_library.txt

# run patch-based segmentation with 4 classes: BG A B C
$bindir/itk_patch_segmentation  --train $rundir/patch_seg_library.txt  --threshold 0.0 --discrete 4 --prob $rundir/patch_seg_prob $rundir/patch_seg_prob_labs.mnc

minccalc -express 'A[0]-A[1]' $rundir/test_sphere_noise_itk_nlm.mnc $rundir/test_sphere.mnc $rundir/test_sphere_noise_itk_anlm_diff.mnc -clob

mean=$(mincstats -q -mean $rundir/test_sphere_noise_itk_anlm_diff.mnc)
var=$(mincstats -q -stddev $rundir/test_sphere_noise_itk_anlm_diff.mnc)

check=$(bc -l <<END
$mean>-$threshold && $mean<$threshold && $var>-$threshold && $var<$threshold
END
)

echo Tests: $mean $var $check

exit $((1-$check))
