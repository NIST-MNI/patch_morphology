#! /bin/sh
set -e
set -x

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

echo Running patch-based grading test
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
             $rundir/grad_test_perfect_sphere_1.mnc

param2xfm -translation 30 0 0 $rundir/grad_shift_x30.xfm -clob

mincresample -nearest -transform $rundir/grad_shift_x30.xfm -use_input $rundir/grad_test_perfect_sphere_1.mnc $rundir/grad_test_perfect_sphere_2.mnc -clob
mincresample -nearest -transform $rundir/grad_shift_x30.xfm -use_input $rundir/grad_test_perfect_sphere_2.mnc $rundir/grad_test_perfect_sphere_3.mnc -clob

#generate training "library"
rm -f $rundir/patch_grad_library.txt
for l in 1 2 3;do
minccalc -express "A[0]>1?$l:0" -byte -labels $rundir/grad_test_perfect_sphere_${l}.mnc $rundir/grad_test_perfect_sphere_${l}_lab.mnc -clob
echo grad_test_perfect_sphere_${l}.mnc,grad_test_perfect_sphere_${l}_lab.mnc >>  $rundir/patch_grade_library.txt
done

random_volume --gauss 1 $rundir/test_perfect_sphere_1.mnc $rundir/test_noise_seg.mnc  --clob

# generate segmentation sample with some noise
minccalc -express 'A[0]+A[1]+A[2]+A[3]' $rundir/test_perfect_sphere_1.mnc $rundir/test_perfect_sphere_2.mnc $rundir/test_perfect_sphere_3.mnc $rundir/test_noise_seg.mnc $rundir/test_seg_sample.mnc -clob

# run patch-based segmentation with 4 classes: BG A B C
$bindir/itk_patch_grading  \
    --exp \
    --train $rundir/patch_seg_library.txt  \
    --patch 1 --search 1 --threshold 0.0 --discrete 4 \
    --prob $rundir/patch_seg_prob_new \
    $rundir/test_seg_sample.mnc \
    $rundir/patch_seg_labs_new.mnc \
    --clob --verbose 

$bindir/itk_patch_morphology \
    --train $rundir/patch_seg_library.txt  \
    --patch 1 --search 1 --threshold 0.0 \
    --discrete 4 \
    --prob $rundir/patch_seg_prob \
    $rundir/test_seg_sample.mnc \
    $rundir/patch_seg_labs.mnc \
    --clob  --verbose 

$bindir/itk_split_labels $rundir/patch_seg_labs_new.mnc \
 $rundir/patch_seg_labs_new_split_%d.mnc --byte 

$bindir/itk_split_labels $rundir/patch_seg_labs.mnc \
 $rundir/patch_seg_labs_split_%d.mnc --byte 


for l in 1 2 3;do
    minccalc -label -express "A[0]-A[1]*$l" $rundir/test_perfect_sphere_${l}_lab.mnc $rundir/patch_seg_labs_new_split_$l.mnc $rundir/patch_seg_labs_new_diff_$l.mnc -clob
    minccalc -label -express "A[0]-A[1]*$l" $rundir/test_perfect_sphere_${l}_lab.mnc $rundir/patch_seg_labs_split_$l.mnc $rundir/patch_seg_labs_diff_$l.mnc -clob
    
    v=$(mincstats -q -sum $rundir/patch_seg_labs_new_diff_$l.mnc)
    if [ "$v" != "0" ];then
        echo "There is a difference in the output , check $rundir/test_perfect_sphere_${l}_lab.mnc $rundir/patch_seg_labs_new_split_$l.mnc"
        exit 1
    fi
    
    v=$(mincstats -q -sum $rundir/patch_seg_labs_diff_$l.mnc)
    if [ "$v" != "0" ];then
        echo "There is a difference in the output , check $rundir/test_perfect_sphere_${l}_lab.mnc $rundir/patch_seg_labs_split_$l.mnc"
        exit 1
    fi
done
