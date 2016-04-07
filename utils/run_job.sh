#!/bin/bash
#
# use the python script run_job.py for more flexibility
# and control this is for legacy machines.

##- options ----------------------------------------------------------
HELP=0
HOSTNAME=`hostname -s|sed s/[[:digit:]]*//g`
JOBNAME=
INITDIR=${SCRATCH}/ves3d
OPTFILE=
WALLTIME=30
NNODE=1
QUEUE=0
RUN=0

##- parsing commandline ----------------------------------------------
OPTIND=1 # reset in case getopts has been used previously in the shell.

while getopts "hH:j:d:o:w:n:qr" opt; do
    case "${opt}" in
	h) HELP=1
           ;;
	H) HOSTNAME=${OPTARG}
           ;;
	j) JOBNAME=${OPTARG}
           ;;
	d) INITDIR=${OPTARG}
           ;;
	o) OPTFILE=${OPTARG}
           ;;
       	w) WALLTIME=${OPTARG}
           ;;
       	n) NNODE=${OPTARG}
           ;;
       	q) QUEUE=1
           ;;
       	r) RUN=1
           ;;
    esac
done

if [ ${HELP} -eq 1 ]; then
    echo "Minimal script to prepare a job, use run_job.py"
    echo "for more flexibility."
    echo
    echo "    -h          : print this message"
    echo "    -H HOST     "
    echo "    -j JOBNAME  "
    echo "    -d INITDIR  : Execuation base directory"
    echo "    -o OPTFILE  : ves3d option file"
    echo "    -w WALLTIME (min)"
    echo "    -n NNODE    "
    echo "    -q QUEUE    "
    echo "    -r RUN      "
    echo
    exit
fi

#- time stamp
STAMP=`date +%m%d.%H%M%S`

echo "Preparing job directory ${INITDIR}"
OPTFILEBASE=$(basename $OPTFILE)
if [ -z ${JOBNAME} ]; then
    JOBNAME=$OPTFILEBASE.${STAMP}
fi
EXECDIR=$INITDIR/$JOBNAME/

# echo "Preparing job file"
if [ -f config/${HOSTNAME}.rc ]; then
    CPCFG="cp config/${HOSTNAME}.rc ${EXECDIR}"
    SRCCFG="source ${HOSTNAME}.rc"
fi

set -x
mkdir -p ${EXECDIR}
cp bin/ves3d ${EXECDIR}
cp ${OPTFILE} $EXECDIR/
${CPCFG}
ln -s ${VES3D_DIR}/precomputed ${EXECDIR}/precomputed
set +x

echo "Preparing job file ${EXECDIR}/${JOBNAME}"
let "WTH=${WALLTIME}/60"
let "WTM=${WALLTIME}%60"
WALLTIME=$WTH:$WTM:00

CMD="mpiexec -np ${NNODE} ves3d -f ${OPTFILEBASE}"

/bin/cat>${EXECDIR}/${JOBNAME} <<EOF
#! /bin/bash

#PBS -l walltime=${WALLTIME}
#PBS -l nodes=${NNODE}
#PBS -m n
#PBS -o localhost:${PBS_W_DIR}
#PBS -j oe

module purge
${SRCCFG}
export VES3D_DIR=${EXECDIR}

${CMD}

EOF

if [ ${QUEUE} -eq 1 ]; then
    cd ${EXECDIR}
    qsub ${JOBNAME}
fi

if [ ${RUN} -eq 1 ]; then
    echo "Running ${CMD}"
    cd ${EXECDIR}
    chmod u+x ${JOBNAME}
    nohup ${CMD} 2>&1 1>${JOBNAME}.out &
fi
