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

# force working with V2 minc files
export MINC_FORCE_V2=1

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


for l in 1 2 3;do
minccalc -express "A[0]>1?$l:0" -byte -labels $rundir/grad_test_perfect_sphere_${l}.mnc $rundir/grad_test_perfect_sphere_${l}_lab.mnc -clob
done



#generate training "library"
rm -f $rundir/patch_grade_library.txt
for l in 1 2 3;do
random_volume --gauss 1 $rundir/grad_test_perfect_sphere_1.mnc $rundir/grad_test_noise_.mnc  --clob
minccalc -express 'A[0]+A[1]+A[2]+A[3]' \
    $rundir/grad_test_noise_.mnc \
    $rundir/grad_test_perfect_sphere_1.mnc \
    $rundir/grad_test_perfect_sphere_2.mnc \
    $rundir/grad_test_perfect_sphere_3.mnc \
    $rundir/grad_sample_${l}.mnc -clob
    
minccalc -labels -express 'A[0]+A[1]+A[2]' \
    $rundir/grad_test_perfect_sphere_1_lab.mnc \
    $rundir/grad_test_perfect_sphere_2_lab.mnc \
    $rundir/grad_test_perfect_sphere_3_lab.mnc \
    $rundir/grad_sample_${l}_seg.mnc  -clob

echo grad_sample_${l}.mnc,grad_sample_${l}_seg.mnc,0,1 >>  $rundir/patch_grade_library.txt
done

# make a dataset with different group
random_volume --gauss 1 $rundir/grad_test_perfect_sphere_1.mnc $rundir/grad_test_noise_.mnc  --clob
minccalc -express 'A[0]+A[1]*1.5+A[2]*1.5+A[3]*1.5' \
    $rundir/grad_test_noise_.mnc \
    $rundir/grad_test_perfect_sphere_1.mnc \
    $rundir/grad_test_perfect_sphere_2.mnc \
    $rundir/grad_test_perfect_sphere_3.mnc \
    $rundir/grad_sample_4.mnc -clob
    
minccalc -labels -express 'A[0]+A[1]+A[2]' \
    $rundir/grad_test_perfect_sphere_1_lab.mnc \
    $rundir/grad_test_perfect_sphere_2_lab.mnc \
    $rundir/grad_test_perfect_sphere_3_lab.mnc \
    $rundir/grad_sample_4_seg.mnc -clob

echo grad_sample_4.mnc,grad_sample_4_seg.mnc,1,2 >>  $rundir/patch_grade_library.txt


# generate noise for new sample
random_volume --gauss 1 $rundir/grad_test_perfect_sphere_1.mnc $rundir/grad_test_noise_seg.mnc  --clob

# generate segmentation sample with some noise
# first sphere is from group  1 
minccalc -express 'A[0]*1.5+A[1]+A[2]+A[3]' \
    $rundir/grad_test_perfect_sphere_1.mnc \
    $rundir/grad_test_perfect_sphere_2.mnc \
    $rundir/grad_test_perfect_sphere_3.mnc \
    $rundir/grad_test_noise_seg.mnc \
    $rundir/test_grad_sample.mnc -clob

# TODO: add --groups
# run patch-based segmentation with 4 classes: BG A B C
$bindir/itk_patch_grading  \
    --exp \
    --train $rundir/patch_grade_library.txt  \
    --patch 1 \
    --search 1 \
    --threshold 0.0 \
    --discrete  4 \
    --grading $rundir/patch_grad_new.mnc \
    $rundir/test_grad_sample.mnc \
    $rundir/patch_grad_labs_new.mnc \
    --clob --verbose 

$bindir/itk_patch_morphology \
    --train $rundir/patch_grade_library.txt  \
    --patch 1 \
    --search 1 \
    --threshold 0.0 \
    --discrete  4 \
    --grading $rundir/patch_grad.mnc \
    $rundir/test_grad_sample.mnc \
    $rundir/patch_grad_labs.mnc \
    --clob  --verbose 

# make sure we properly found all the labels and calculated gradings 
itk_label_stats $rundir/patch_grad_labs_new.mnc --volume $rundir/patch_grad_new.mnc > $rundir/grading_label_stats_new.csv
itk_label_stats $rundir/patch_grad_labs.mnc --volume $rundir/patch_grad.mnc > $rundir/grading_label_stats.csv

# ground truth for volumes of labels:
cat - >$rundir/grading_truth.csv <<END
id,volume,mx,my,mz,val
1,936,-30,0,0,1
2,936,0,0,0,0
3,936,30,0,0,0
END

if ! diff -q $rundir/grading_label_stats_new.csv $rundir/grading_truth.csv;then
echo "Unexpected values in $rundir/grading_label_stats_new.csv"
exit 1
fi

if ! diff -q $rundir/grading_label_stats.csv $rundir/grading_truth.csv;then
echo "Unexpected values in $rundir/grading_label_stats.csv"
exit 1
fi

