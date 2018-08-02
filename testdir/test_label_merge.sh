#! /bin/bash
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
             $rundir/test_perfect_sphere_1.mnc

param2xfm -translation 30 0 0 $rundir/shift_x30.xfm -clob

mincresample -nearest -transform $rundir/shift_x30.xfm -use_input $rundir/test_perfect_sphere_1.mnc $rundir/test_perfect_sphere_2.mnc -clob
mincresample -nearest -transform $rundir/shift_x30.xfm -use_input $rundir/test_perfect_sphere_2.mnc $rundir/test_perfect_sphere_3.mnc -clob

for l in 1 2 3;do
    minccalc -express "A[0]>1?$l:0" -byte -label $rundir/test_perfect_sphere_$l.mnc $rundir/test_perfect_sphere_lab_$l.mnc -clob
done
minccalc -float -express '0.01' $rundir/test_perfect_sphere_1.mnc $rundir/test_bg.mnc -clob

# generate segmentation with multiple labels
minccalc -byte -label -express 'A[0]+A[1]+A[2]' \
    $rundir/test_perfect_sphere_lab_1.mnc \
    $rundir/test_perfect_sphere_lab_2.mnc \
    $rundir/test_perfect_sphere_lab_3.mnc \
    $rundir/test_perfect_sphere_labs.mnc -clob

# do the same with merge_labels on non-binary volume
$bindir/itk_merge_labels \
  0 $rundir/test_bg.mnc  \
  1 $rundir/test_perfect_sphere_1.mnc  \
  2 $rundir/test_perfect_sphere_2.mnc  \
  3 $rundir/test_perfect_sphere_3.mnc  \
    $rundir/test_perfect_sphere_labs_merged.mnc --clob --byte


    
# make sure volumes are the same
minccalc -byte -label -express 'abs(A[0]-A[1])>0.1?1:0' $rundir/test_perfect_sphere_labs.mnc $rundir/test_perfect_sphere_labs_merged.mnc \
 $rundir/test_perfect_sphere_labs_merged_diff.mnc -clob

v=$(mincstats -sum -q $rundir/test_perfect_sphere_labs_merged_diff.mnc)

if [ "$v" != "0" ];then
    echo "Unexpected results after merging, check $rundir/test_perfect_sphere_labs_merged.mnc"
    exit 1
fi
 
