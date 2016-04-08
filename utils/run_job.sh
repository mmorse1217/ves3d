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

OPTFILEBASE=$(basename $OPTFILE)
if [ -z ${JOBNAME} ]; then
    JOBNAME=$OPTFILEBASE.${STAMP}
fi
EXECDIR=$INITDIR/$JOBNAME/
echo "Preparing job directory ${EXECDIR}"

if [ -f config/${HOSTNAME}.rc ]; then
    CPCFG="cp config/${HOSTNAME}.rc ${EXECDIR}/"
    SRCCFG="source ${HOSTNAME}.rc"
fi

set -x
mkdir -p ${EXECDIR}
cp bin/ves3d ${EXECDIR}/
cp ${OPTFILE} ${EXECDIR}/
${CPCFG}
cp utils/epilogue.sh ${EXECDIR}/
ln -s ${VES3D_DIR}/precomputed ${EXECDIR}/precomputed
set +x

echo "Preparing job file ${EXECDIR}/${JOBNAME}"
let "WTH=${WALLTIME}/60"
let "WTM=${WALLTIME}%60"
WALLTIME=$WTH:$WTM:00

CMD="mpiexec --bind-to-core --bynode -np \${PBS_NUM_NODES}"
CMD="${CMD} -x PATH -x LD_LIBRARY_PATH ves3d -f ${OPTFILEBASE} 2>&1 1>${JOBNAME}.out"

/bin/cat>${EXECDIR}/${JOBNAME} <<EOF
#!/bin/bash

#PBS -l walltime=${WALLTIME}
#PBS -l nodes=${NNODE}
#PBS -l epilogue=epilogue.sh
#PBS -m n
#PBS -o localhost:\${PBS_O_WORKDIR}
#PBS -j oe

echo ----------------------------------------------------------------------
echo PBS: qsub is running on \${PBS_O_HOST}
echo PBS: originating queue is \${PBS_O_QUEUE}
echo PBS: executing queue is \${PBS_QUEUE}
echo PBS: submission directory is \${PBS_O_WORKDIR}
echo PBS: execution mode is \${PBS_ENVIRONMENT}
echo PBS: job identifier is \${PBS_JOBID}
echo PBS: job name is \${PBS_JOBNAME}
echo PBS: node file is \${PBS_NODEFILE}
echo PBS: current home directory is \${PBS_O_HOME}
echo PBS: PATH = \${PBS_O_PATH}
echo ----------------------------------------------------------------------
echo PBS: NP = \${PBS_NP}
echo PBS: NUM_PPN = \${PBS_NUM_PPN}
echo ----------------------------------------------------------------------

cd \${PBS_O_WORKDIR}
source ${MODULESHOME}/init/\$(basename \$SHELL)

${SRCCFG}
export VES3D_DIR=${EXECDIR}
export OMP_NUM_THREADS=\${PBS_NUM_PPN}
env
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
