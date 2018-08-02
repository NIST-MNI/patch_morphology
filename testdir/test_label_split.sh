#! /bin/bash
set -e 

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
    minccalc -labels -express "A[0]>0.5?$l:0" -byte $rundir/test_perfect_sphere_$l.mnc $rundir/test_perfect_sphere_lab_$l.mnc -clob
done

# generate segmentation with multiple labels
minccalc -byte -labels -express 'A[0]+A[1]+A[2]' \
    $rundir/test_perfect_sphere_lab_1.mnc \
    $rundir/test_perfect_sphere_lab_2.mnc \
    $rundir/test_perfect_sphere_lab_3.mnc \
    $rundir/test_perfect_sphere_labs.mnc -clob

# try with simple threshold
$bindir/itk_split_labels \
 --byte --clob \
    $rundir/test_perfect_sphere_labs.mnc  \
    $rundir/test_perfect_sphere_labs_split_hard_%d.mnc

# compare volumes
for l in 1 2 3;do
    v1=$(mincstats -q -floor 0.5 -vol $rundir/test_perfect_sphere_lab_$l.mnc)
    v2=$(mincstats -q -floor 0.5 -vol $rundir/test_perfect_sphere_labs_split_hard_$l.mnc)
    
    if [ "$v1" != "$v2" ];then
        echo "Volumes $rundir/test_perfect_sphere_lab_$l.mnc and $rundir/test_perfect_sphere_labs_split_hard_$l.mnc mismatch"
        exit 1
    fi
done
    
# try with antialias
$bindir/itk_split_labels \
 --float --aa --clob \
    $rundir/test_perfect_sphere_labs.mnc  \
    $rundir/test_perfect_sphere_labs_split_soft_%d.mnc

# compare volumes
for l in 1 2 3;do
    v1=$(mincstats -q -floor 0.5 -vol $rundir/test_perfect_sphere_lab_$l.mnc)
    v2=$(mincstats -q -floor 0.0 -vol $rundir/test_perfect_sphere_labs_split_soft_$l.mnc)
    
    if [ "$v1" != "$v2" ];then
        echo "Volumes $rundir/test_perfect_sphere_lab_$l.mnc and $rundir/test_perfect_sphere_labs_split_soft_$l.mnc mismatch"
        exit 1
    fi
done

# try with antialias and expit function
$bindir/itk_split_labels \
 --float --aa --clob --expit 0.5 \
    $rundir/test_perfect_sphere_labs.mnc  \
    $rundir/test_perfect_sphere_labs_split_soft_expit_%d.mnc

# compare volumes
for l in 1 2 3;do
    v1=$(mincstats -q -floor 0.5 -vol $rundir/test_perfect_sphere_lab_$l.mnc)
    v2=$(mincstats -q -floor 0.5 -vol $rundir/test_perfect_sphere_labs_split_soft_expit_$l.mnc)
    
    if [ "$v1" != "$v2" ];then
        echo "Volumes $rundir/test_perfect_sphere_lab_$l.mnc and $rundir/test_perfect_sphere_labs_split_soft_expit_$l.mnc mismatch"
        exit 1
    fi
done

# try with missing label 
$bindir/itk_split_labels \
 --byte --clob --missing 5 \
    $rundir/test_perfect_sphere_labs.mnc  \
    $rundir/test_perfect_sphere_labs_split_missing_%d.mnc

# should have zero volume
v2=$(mincstats -q -floor 0.5 -vol $rundir/test_perfect_sphere_labs_split_missing_4.mnc)    
if [ "0" != "$v2" ];then
    echo "Volume $rundir/test_perfect_sphere_labs_split_missing_4.mnc is not empty"
    exit 1
fi
