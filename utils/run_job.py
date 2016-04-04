#!/usr/bin/env python
'''Manages the execution of a job. This creates pbs job file, job-specific
option file, etc.'''

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
import sys
import time
import yaml
import subprocess as sp

class CustomArgFormatter(ap.ArgumentDefaultsHelpFormatter,
                         ap.RawDescriptionHelpFormatter):
    pass

class Job(object):
    def __init__(self,**kwargs):
        self._opts = self._parse_cl(**kwargs)

        if self.host is None:
            host = sp.Popen(['hostname', '-s'], stdout=sp.PIPE).communicate()
            self.host = host[0].strip()
        self._host_opts = self._load_host_opts()

        self.set_defaults()

        self.init_basedir = os.path.expandvars(self.init_basedir)
        self.init_dir     = os.path.join(self.init_basedir,self.job_name+'/')

        self.exec_name    = os.path.abspath(self.exec_name)
        self.optfile      = os.path.abspath(self.optfile)

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

        clp.add_argument('--host', '-H', help='Hostname')
        clp.add_argument('--hostfile', '-F', help='Hostname')

        #queue
        clp.add_argument('--account'   , '-A' , help='Account number')
        clp.add_argument('--job-name'  , '-N' , help='Job name')
        clp.add_argument('--queue'     , '-q' , help='Destination queue')

        #logging/notification
        clp.add_argument('--outfile'   , '-o' , help='Output file',default='localhost:${PBS_O_WORKDIR}/')
        clp.add_argument('--errfile'   , '-e' , help='Error file')
        clp.add_argument('--join'      , '-j' , help='Join out/err files' ,default='oe')
        clp.add_argument('--mail'      ,        help='Mail options [begin, abort, end]', default='n', choices=('b','e','a','s','n'))
        clp.add_argument('--userlist'  ,        help='List of users (to notify)')

        #resources
        clp.add_argument('--nodes'     , '-n' , help='Number of nodes')
        clp.add_argument('--ppn'       , '-p' , help='Processor per node')
        clp.add_argument('--walltime'  , '-w' , help='Wall clock time (use suffix h for hour, otherwise interpreted as minutes)')
        clp.add_argument('--memory'    , '-m' , help='Memory per node (GB)')
        clp.add_argument('--resources' , '-l' , help='Extra resources (comma separated)')

        #setup (not for pbs)
        clp.add_argument('--no-bin'       , help='Skip copying binary file and symlink', action='store_true')
        clp.add_argument('--no-manifest'  , help='Skip manifest file', action='store_true')

        #executable
        clp.add_argument('--exec-name'    , '-E' , help='Executable name (required)',required=True)
        clp.add_argument('--optfile'      , '-I' , help='Input option file for the executable)')
        clp.add_argument('--init-basedir' , '-i' , help='Init directory')
        clp.add_argument('--vars'         , '-v' , help='Variable list (comma separated)')
        clp.add_argument('--modules'      ,        help='modules to load')

        if callback is not None: callback(clp)
        return vars(clp.parse_args())

    def _load_host_opts(self):
        hostfile = self.hostfile
        if hostfile is None:
            hostfile = 'config/%s.yml' % self.host

        with open(hostfile, 'r') as stream:
            try:
                mopts = yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)

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
                self.job_name = '%s.%s' % (jname,stamp)

        if self.walltime is not None:
            h = m = 0
            if self.walltime[-1]=='h':
                h = int(self.walltime[:-1])
            else:
                m = int(self.walltime)
            self.walltime = '%02d:%02d:00' % (h,m)

        if self.memory is not None:
            self.memory += 'G'

        if self.ppn is not None:
            self.nodes +=':ppn=%s'%self.ppn

    def pbs_args(self):
        args = []
        def append_args(lst):
            for k1,k2 in lst:
                args.append((k2,getattr(self,k1)))

        queue     = [('account','-A '), ('job_name','-N '), ('queue','-q ')]
        append_args(queue)

        logging   = [('outfile','-o '), ('errfile','-e '), ('join','-j '),
                     ('mail','-m '), ('userlist','-M ')]
        append_args(logging)

        resources = [('nodes','-l node='), ('walltime','-l walltime='),
                     ('memory','-l mem='), ('resources','-l ')]
        append_args(resources)

        exe       = [('init_dir','-d '), ('vars','-v ')]
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

    def pbs_job_header(self):

        pbs  = ['#!/bin/bash\n']
        args = self.pbs_args()

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

        cmd  = '${MPI_ROOT}/bin/mpiexec -np ${NT} -pernode'
        cmd += ' -x PATH -x LD_LIBRARY_PATH -x OMP_NUM_THREADS'
        cmd += ' -wdir ${PBS_O_WORKDIR} ${PBS_O_WORKDIR}/%s -f %s' % (efile,ofile)

        mods = 'module purge\n'
        for m in self.modules:
            mods += 'module load %s\n' % m

        if self.host=='octane':
            cmds = [
                'cd ${PBS_O_WORKDIR}',
                '',
                '#set number of OpenMP threads per node, NT=number of MPI nodes (not threads)',
                'export OMP_NUM_THREADS=${PBS_NUM_PPN}',
                'export NT=$((PBS_NP/PBS_NUM_PPN))',
                '',
                mods,
                'CMD='+cmd,
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
        content  = self.pbs_job_header()
        content += ['\n']
        content += self.pbs_log_header()
        content += ['\n']
        content += self.exec_cmd(fnames)
        content += ['\n','#EOF']

        with open(jname,'w') as fh:
            for l in content:
                fh.write(l+'\n')

        return jname

    def prepare_dir(self):
        print('preparing init ldir', self.init_dir)
        os.makedirs(self.init_dir,0770)

        #executable
        execdir  = os.path.dirname(self.exec_name)
        execname = os.path.basename(self.exec_name)
        execname = os.path.join(self.init_dir,execname)

        if self.no_bin:
            os.symlink(self.exec_name,execname)
        else:
            sp.call(['cp',self.exec_name,execname])

        #option file
        optfile = os.path.basename(self.optfile)
        optfile = os.path.join(self.init_dir,optfile)
        sp.call(['cp',self.optfile,optfile])

        #precomputed
        precomp = os.path.expandvars('${VES3D_DIR}')
        precomp = os.path.join(precomp,'precomputed')
        assert os.path.isdir(precomp), "precomputed file '%s' does not exist" % precomp
        os.symlink(precomp,os.path.join(self.init_dir,'precomputed'))

        #file manifest
        if not self.no_manifest:
            mfile = os.path.join(self.init_dir,'manifest')
            with open(mfile, 'w') as fh:
                sp.Popen('hg summary',shell=True,stdout=fh,stderr=fh)
                sp.Popen('hg status',shell=True,stdout=fh,stderr=fh)

        return dict(execname=execname,optfile=optfile)

    def prepare(self):
        fnames = self.prepare_dir()
        jname  = self.prepare_job_file(fnames)

        return fnames, jname
    def submit(self):
        fnames, jname = self.prepare()

        print('submitting %s (%s)' % (self.job_name, jname))
        os.chdir(self.init_dir)
        sp.call(['qsub', jname])

Job().submit()
