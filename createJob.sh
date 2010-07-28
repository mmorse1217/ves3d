#!/bin/sh
#input EXECUTABLENAME [JOBFILE] [OMP_NUM_THREADS] [LABEL]

#--- general options ------------------------------------
USEREMAIL="rahimian@gatech.edu"
GRANTNUMBER=""

MACHINE=$(hostname)
WALLTIME=24:00:00
CODEVERSION="VES3D_$(git tag)"
SRCDIR=$VES3D_DIR
JOBFILE=$2

#- openmp
NP=1
OPENMP=y
NUMTHREADS=$3

if [ -z $NUMTHREADS ]; then
    NUMTHREADS=1
fi

if [ -n $OPENMP ]; then
    OMPTHREADS="export OMP_NUM_THREADS=$NUMTHREADS"
    OPENMPLABEL="omp$NUMTHREADS"
    NPREQUESTED=$NP
fi

#- label
USETIMELABEL=y
LABEL=$4

MON=`date +%b`
DAY=`date +%d`
HRM=`date +%H%M`

if [ $LABEL ]; then
    LABEL="$LABEL."
fi

if [ $USETIMELABEL ]; then
    LABEL="$MON$DAY.$HRM.$LABEL"
fi
LABEL="$LABEL$OPENMPLABEL.p$NPREQUESTED"

#- filename options
EXECUTABLE=$1
EXECPATH=`dirname $EXECUTABLE`
EXECUTABLE=${EXECUTABLE##*/}
EXECUTABLE=${EXECUTABLE%.*}

SCRATCHDIR=$SRCDIR/scratch/

TARGETEXECUTABLE="$EXECUTABLE.$LABEL.exe"
OUTFILE="$EXECUTABLE.$LABEL.out"

#- vesicles options
if [ $JOBFILE ]; then
    JOBFILE=$JOBFILE.$EXECUTABLE.$LABEL
    exec > $JOBFILE
fi

#--- generating the job ---------------------------------
echo "#PBS -M $USEREMAIL" 
echo "#PBS -m ae"
if [ $GRANTNUMBER ]; then
    echo "#PBS -A $GRANTNUMBER"
fi
echo "#PBS -l walltime=$WALLTIME"
echo "#PBS -l nodes=$NPREQUESTED"
echo "#PBS -N $TARGETEXECUTABLE"
echo
echo "cp $SRCDIR/$EXECPATH/$EXECUTABLE.exe $SCRATCHDIR/$TARGETEXECUTABLE"
echo "$OMPTHREADS"
echo 
echo "cd $SCRATCHDIR"
echo "./$TARGETEXECUTABLE > $OUTFILE"
echo "mv $OUTFILE $SRCDIR/results"
