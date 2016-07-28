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

echo Running nlm test
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
make_phantom -ellipse \
             -continuous -short \
             -start -50 -50 -50 \
             -nelements 50 50 50 \
             -background 0 \
             -step 2 2 2 \
             -fill_value 100 \
             -edge_value 50  \
             -width 20 20 20 \
             -clob \
             -center 0 0 0  \
             $rundir/test_sphere.mnc
             
random_volume --gauss 10 $rundir/test_sphere.mnc $rundir/test_noise.mnc  --clob
minccalc -express 'A[0]+A[1]' $rundir/test_sphere.mnc $rundir/test_noise.mnc $rundir/test_sphere_noise.mnc -clob
$bindir/itk_minc_nonlocal_filter --anlm $rundir/test_sphere_noise.mnc $rundir/test_sphere_noise_itk_anlm.mnc  --flat --clob 
minccalc -express 'A[0]-A[1]' $rundir/test_sphere_noise_itk_nlm.mnc $rundir/test_sphere.mnc $rundir/test_sphere_noise_itk_anlm_diff.mnc -clob

mean=$(mincstats -q -mean $rundir/test_sphere_noise_itk_anlm_diff.mnc)
var=$(mincstats -q -stddev $rundir/test_sphere_noise_itk_anlm_diff.mnc)

check=$(bc -l <<END
$mean>-$threshold && $mean<$threshold && $var>-$threshold && $var<$threshold
END
)

echo Tests: $mean $var $check

exit $((1-$check))
