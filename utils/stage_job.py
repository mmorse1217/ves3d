#!/usr/bin/env python
'''Manages the execution of a job. This creates pbs job file, job-specific
option file, etc.

This reads a hostfile to set the default values based on the machine.
The machine dependent options can be set there.
'''

from __future__ import absolute_import, division, print_function

__author__    = 'Abtin Rahimian'
__email__     = 'arahimian@acm.org'
__status__    = 'prototype'
__revision__  = '$Revision$'
__date__      = '$Date$'
__tags__      = '$Tags$'
__copyright__ = 'Copyright (c) 2016, Abtin Rahimian'

import argparse as ap
import datetime
import os
import re
import sys
import time
import subprocess as sp

class CustomArgFormatter(ap.ArgumentDefaultsHelpFormatter,
                         ap.RawDescriptionHelpFormatter):
    pass

class Job(object):
    def __init__(self,**kwargs):
        self._opts = self._parse_cl(**kwargs)

        self.host = os.path.expandvars(self.host)
        self._host_opts = self._load_host_opts()

        self.set_defaults()

        self.init_basedir = os.path.expandvars(self.init_basedir)
        self.init_dir     = os.path.join(self.init_basedir,self.job_name+'/')

        self.exec_name    = os.path.abspath(self.exec_name)
        self.optfile      = os.path.abspath(self.optfile)
        self.epilogue     = os.path.abspath(self.epilogue)

    def __getattr__(self,name):
        return self._opts[name]

    def _parse_cl(self,callback=None,predoc=None,postdoc=None):
        # commandline parser
        clp = ap.ArgumentParser(
            formatter_class=CustomArgFormatter,
            description=predoc,epilog=postdoc,
            )

        clp.add_argument('--queue-manager', '-Q',
                         help='Queue management system', choices=('svg','torque'),
                         default='torque')

        clp.add_argument('--host'    , '-H', help='Hostname (if an evironment variable, it is expanded)', default='${VES3D_PLATFORM}')
        clp.add_argument('--hostfile', '-F', help='Hostfile to read machine defaults (if unset, it is assumed to be <host>.yml)')

        #queue
        clp.add_argument('--account'   , '-A' , help='Account number')
        clp.add_argument('--job-name'  , '-N' , help='Job name (default to optfile + timestamp)')
        clp.add_argument('--queue'     , '-q' , help='Destination queue')
        clp.add_argument('--no-stage'  ,        help='Skip submitting the job', action='store_true')

        #logging/notification
        clp.add_argument('--outfile'   , '-o' , help='Output file', default='localhost:${PBS_O_WORKDIR}/')
        clp.add_argument('--errfile'   , '-e' , help='Error file (by default joined to outfile)')
        clp.add_argument('--join'      , '-j' , help='Join out/err files', default='oe')
        clp.add_argument('--mail'      ,        help='Mail options [begin, end, abort, none]', default='ea', choices=('b','e','a','s','n'))
        clp.add_argument('--userlist'  ,        help='List of users to notify')

        #resources
        clp.add_argument('--nodes'     , '-n' , help='Number of nodes',type=int)
        clp.add_argument('--cpu'       , '-c' , help='Number of cpus per node',type=int)
        clp.add_argument('--ppn'       , '-p' , help='Number of mpi processes per node',type=int,default=1)
        clp.add_argument('--threads'   , '-t' , help='Number of threads per process',type=int)
        clp.add_argument('--walltime'  , '-w' , help='Wall clock time (use suffix h for hour, otherwise interpreted as minutes)')
        clp.add_argument('--memory'    , '-m' , help='Memory per node (GB if no unit)')
        clp.add_argument('--resources' , '-l' , help='Extra resources (copied verbatim)')
        clp.add_argument('--epilogue'  ,        help='Epilogue script to run after code', default='utils/epilogue.sh')

        #setup (not for pbs)
        clp.add_argument('--no-bin'       , help='Skip copying binary file and symlink', action='store_true')
        clp.add_argument('--no-manifest'  , help='Skip manifest file', action='store_true')

        #executable
        clp.add_argument('--exec-name'    , '-E' , help='Executable name', default='bin/ves3d')
        clp.add_argument('--optfile'      , '-I' , help='Input option file for the executable', required=True)
        clp.add_argument('--init-basedir' , '-i' , help='Init directory--the location where the job is setup, submitted, and output files stored (environment variables are expanded)',default='${SCRATCH}/ves3d/')
        clp.add_argument('--vars'         , '-v' , help='Variable list (copied verbatim)')
        clp.add_argument('--modules'      ,        help='modules to load',nargs='*')

        if callback is not None: callback(clp)
        return vars(clp.parse_args())

    def _load_host_opts(self):
        hostfile = self.hostfile
        if hostfile is None:
            hostfile = 'config/%s.yml' % self.host

        if not os.path.exists(hostfile):
            return dict()

        print('loading host file %s' % hostfile)
        fh = open(hostfile, 'r')
        import yaml
        #if missing, use pip install PyYAML
        mopts = yaml.load(fh)
        fh.close()

        for k,v in mopts.items():
            nk = k.replace('-','_')
            if nk!=k:
                mopts[nk]=v
                del mopts[k]

        return mopts

    def set_defaults(self):
        stamp = datetime.datetime.fromtimestamp(time.time())
        stamp = stamp.strftime('%m%d.%H%M%S')

        for k,v in self._host_opts.items():
            if getattr(self,k) is None:
                setattr(self,k,v)

        if self.job_name is None and self.optfile is not None:
            jname = os.path.basename(self.optfile)
            if len(jname):
                self.job_name = '%s.%s.pbs' % (jname,stamp)

        if self.walltime is not None:
            h = m = 0
            if self.walltime[-1]=='h':
                h = int(self.walltime[:-1])
            else:
                m = int(self.walltime)
            self.walltime = '%02d:%02d:00' % (h,m)

        if self.memory is not None:
            if re.search('[kKmMgG]', self.memory) is None:
                self.memory += 'G'

        if self.cpu is not None:
            self.nodes ='%d:ppn=%d'% (self.nodes,self.cpu)

        if len(self.epilogue)==0: self.epilogue=None

    def pbs_args(self,files):
        args = []
        def append_args(lst):
            for k1,k2 in lst:
                args.append((k2,getattr(self,k1)))

        queue     = [('account','-A '), ('job_name','-N '), ('queue','-q ')]
        append_args(queue)

        logging   = [('outfile','-o '), ('errfile','-e '), ('join','-j '),
                     ('mail','-m '), ('userlist','-M ')]
        append_args(logging)

        resources = [('nodes','-l nodes='), ('walltime','-l walltime='),
                     ('memory','-l mem='), ('resources','-l ')]
        append_args(resources)

        epilogue = files.get('epilogue', None)
        if  epilogue is not None:
            args.append(('-l epilogue=', os.path.basename(epilogue)))

        exe       = [('vars','-v ')]
        append_args(exe)

        return args

    def pbs_log_header(self):
        header = [
            '----------------------------------------------------------------------',
            'PBS: qsub is running on ${PBS_O_HOST}',
            'PBS: originating queue is ${PBS_O_QUEUE}',
            'PBS: executing queue is ${PBS_QUEUE}',
            'PBS: submission directory is ${PBS_O_WORKDIR}',
            'PBS: execution mode is ${PBS_ENVIRONMENT}',
            'PBS: job identifier is ${PBS_JOBID}',
            'PBS: job name is ${PBS_JOBNAME}',
            'PBS: node file is ${PBS_NODEFILE}',
            'PBS: current home directory is ${PBS_O_HOME}',
            'PBS: PATH = ${PBS_O_PATH}',
            '----------------------------------------------------------------------',
            'PBS: NP = ${PBS_NP}',
            'PBS: NUM_PPN = ${PBS_NUM_PPN}',
            '----------------------------------------------------------------------',
            ]

        header = ['echo '+l for l in header]
        return header

    def pbs_job_header(self,files):

        pbs  = ['#!/bin/bash\n']
        args = self.pbs_args(files)

        for k,v in args:
            if k is None: pbs.append('')
            if v is not None: pbs.append('#PBS %s%s'%(k,v))

        return pbs

    def exec_cmd(self,fnames):

        assert os.path.exists(fnames['execname'])
        assert os.path.exists(fnames['optfile'])

        edir  = os.path.dirname(fnames['execname'])
        efile = os.path.basename(fnames['execname'])
        ofile = os.path.basename(fnames['optfile'])
        assert os.path.samefile(edir,self.init_dir), 'inconsistent dir name'

        mods = 'module purge\n'
        for m in self.modules:
            mods += 'module load %s\n' % m
        mods += 'module list\n'

        if self.threads is None:
            if self.cpu is None:
                omp_threads='OMP_NUM_THREADS=$((PBS_NUM_PPN/%d))'%self.ppn
            else:
                omp_threads='OMP_NUM_THREADS=%d'% (self.cpu/self.ppn)
        else:
            omp_threads='OMP_NUM_THREADS=%d'%self.threads

        if self.host=='mercer':
            # bind-to-core to avoid switching CPU of an mpi process
            # need bynode so each mpi process is in a different node (works for openmpi)
            #mpirun --bind-to-core --bynode -np $PBS_NP
            cmd  = 'mpiexec --bind-to-core -bynode -loadbalance -np ${NP}'
            cmd += ' -x OMP_NUM_THREADS -x PATH -x LD_LIBRARY_PATH ./%s -f %s' % (efile,ofile)

            cmds = [
                'cd ${PBS_O_WORKDIR}\n',
                'export VES3D_DIR=%s' % self.init_dir,
                'export %s'%omp_threads,
                'NP=$((PBS_NUM_NODES*%d))'%self.ppn,
                '',
                mods,
                'echo Environment variables:\nenv\n\n',
                'CMD="%s"'%cmd,
                'echo running ${CMD}',
                '${CMD}',
                ]
        else:
            raise 'Do not know how to run executable on %s' % self.host

        return cmds

    def prepare_job_file(self,fnames):

        jname = os.path.join(self.init_dir,self.job_name)
        print('preparing job file %s' % jname)

        content  = []
        content  = self.pbs_job_header(fnames)
        content += ['\n']
        content += self.pbs_log_header()
        content += ['\n']
        content += self.exec_cmd(fnames)
        content += ['\n','#EOF']

        fh = open(jname,'w')
        for l in content: fh.write(l+'\n')
        fh.close()

        return jname

    def prepare_dir(self):
        print('preparing init dir', self.init_dir)
        os.makedirs(self.init_dir,0770)
        files = dict()
        #executable
        execdir  = os.path.dirname(self.exec_name)
        execname = os.path.basename(self.exec_name)
        execname = os.path.join(self.init_dir,execname)
        files['execname'] = execname

        if self.no_bin:
            os.symlink(self.exec_name,execname)
        else:
            sp.call(['cp',self.exec_name,execname])

        #option file
        optfile = os.path.basename(self.optfile)
        optfile = os.path.join(self.init_dir,optfile)
        sp.call(['cp',self.optfile,optfile])
        files['optfile'] = optfile

        #precomputed
        precomp = os.path.expandvars('${VES3D_DIR}')
        precomp = os.path.join(precomp,'precomputed')
        assert os.path.isdir(precomp), "precomputed file '%s' does not exist" % precomp
        os.symlink(precomp,os.path.join(self.init_dir,'precomputed'))

        #file manifest
        if not self.no_manifest:
            mfile = os.path.join(self.init_dir,'manifest')
            fh = open(mfile, 'w')
            sp.Popen(['hg', 'summary'],shell=False,stdout=fh,stderr=fh).wait()
            sp.Popen(['hg', 'status'],shell=False,stdout=fh,stderr=fh).wait()
            fh.close()

        #epilogue file
        if self.epilogue is not None:
            epilogue = os.path.basename(self.epilogue)
            epilogue = os.path.join(self.init_dir,epilogue)
            sp.call(['cp',self.epilogue,epilogue])
            os.chmod(epilogue,0755)
            files['epilogue'] = epilogue
        return files

    def prepare(self):
        fnames = self.prepare_dir()
        jname  = self.prepare_job_file(fnames)

        return fnames, jname

    def stage(self):
        fnames, jname = self.prepare()

        if not self.no_stage:
            print('submitting %s (%s)' % (self.job_name, jname))
            os.chdir(self.init_dir)
            sp.call(['qsub', jname])

Job().stage()
