#!/bin/sh
#input EXECUTABLENAME [JOBFILE_SUFFIX] [OMP_NUM_THREADS] [LABEL_BASE]

#--- general options ------------------------------------
useremail="rahimian@gatech.edu"
grantnumber=""

machine=$(hostname)
walltime=24:00:00
codeversion="VES3D_$(git tag)"
srcdir=$VES3D_DIR
jobfile=$2

#- openmp
np=1
openmp=y
numthreads=$3

if [ -z $numthreads ]; then
    numthreads=1
fi

if [ -n $openmp ]; then
    ompthreads="export OMP_NUM_THREADS=$numthreads"
    openmplabel="omp$numthreads"
    nprequested=$np
fi

#- label
usetimelabel=y
label=$4

mon=`date +%b`
day=`date +%d`
hrm=`date +%H%M`

if [ $label ]; then
    label=$label.
fi

if [ $usetimelabel ]; then
    label=$mon$day.$hrm.$label
fi
label=$label$openmplabel.p$nprequested

#- filename generation
executable=$1
execpath=`dirname $executable`
executable=${executable##*/}
executable=${executable%.*}

scratchdir=$srcdir/scratch/

targetexecutable=$label.$executable.exe
outfile=$label.$executable.out

#- checking the output format
if [ $jobfile ]; then
    jobfile=$label.$executable.$jobfile
    exec > $jobfile  #redirecting stdout to the file
fi

#--- generating the job file ------------------------------
echo "#PBS -M $useremail" 
echo "#PBS -m ae"
if [ $grantnumber ]; then
    echo "#PBS -A $grantnumber"
fi
echo "#PBS -l walltime=$walltime"
echo "#PBS -l nodes=$nprequested"
echo "#PBS -N $targetexecutable"
echo
echo "cp $srcdir/$execpath/$executable.exe $scratchdir/$targetexecutable"
echo "$ompthreads"
echo 
echo "cd $scratchdir"
echo "./$targetexecutable > $outfile"
echo "mv $outfile $srcdir/results"
