#!/bin/bash
#
# commandline option:
#   executablename [jobfile_suffix] [walltime (min)] [OMP_NUM_THREADS] [label]

######################################################################
######################################################################
#                                                                    # 
#    This is the base for the job generating script. Try not to      #         
#    edit this file. Copy and presonalize the last portion of        #
#    this file to another (local) file and source this file          # 
#    there.                                                          #
#                                                                    #
######################################################################
######################################################################

#--- general options -------------------------------------------------
machine=$(hostname)
walltime=$3
: ${walltime:=30} #min
codeversion="VES3D_$(git tag)"
srcdir=$VES3D_DIR
jobfile=$2

#- converting the time to the desired format
let "wth=$walltime/60"
let "wtm=$walltime%60"
walltime=$wth:$wtm:00

#- openmp
np=1
openmp=y
numthreads=$4
: ${numthreads:=1}

if [ -n $openmp ]; then
    ompthreads="export OMP_NUM_THREADS=$numthreads"
    openmplabel="omp$numthreads"
    let "nprequested=np"
fi

#- label
usetimelabel=y
label=$5

mon=`date +%b`
day=`date +%d`
hrm=`date +%H%M`

if [ $label ]; then
    label=$label.
fi

if [ $usetimelabel ]; then
    label=$mon$day.$hrm.$label
fi
label=$label$openmplabel.p$nprequested.$machine

#- filename generation
executable=$1
execpath=`dirname $executable`
executable=${executable##*/}
executable=${executable%.*}

scratchdir=$SCRATCH

targetexecutable=$label.$executable.exe
outfile=$label.$executable.out

#- checking the output format
if [ $jobfile ]; then
    jobfile=$label.$executable.$jobfile
    exec > $jobfile  #redirecting stdout to the file
fi
  
if [ ! -d $scratchdir/precomputed ]; then
    echo `ln -s $VES3D_DIR/precomputed/ $scratchdir/precomputed`
fi

#- Moving the executable to the scratch directory
echo `cp $srcdir/$execpath/$executable.exe $scratchdir/$targetexecutable`

######################################################################
#  Copy and uncomment one of the following
######################################################################
                                                                     
# #!/bin/bash
#
# useremail="rahimian@gatech.edu"
# grantnumber=""
#
# source createJob.sh
#
# echo "#PBS -M $useremail" 
# if [ $grantnumber ]; then
#     echo "#PBS -A $grantnumber"
# fi
# echo "#PBS -l walltime=$walltime"
# echo "#PBS -l nodes=$nprequested"
# echo "#PBS -N $targetexecutable"
# echo "#PBS -m a"
# echo
# echo "$ompthreads"
# echo 
# echo "cd $scratchdir"
# echo "./$targetexecutable > $outfile"
# echo "mv $outfile $srcdir/results/"
#
######################################################################

