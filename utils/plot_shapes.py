#!/usr/bin/env python
'''plot the shapes stored in the input file'''

from __future__ import absolute_import, division, print_function

__author__    = 'Abtin Rahimian'
__email__     = 'arahimian@acm.org'
__status__    = 'prototype'
__revision__  = '$Revision$'
__date__      = '$Date$'
__tags__      = '$Tags$'
__copyright__ = 'Copyright (c) 2015, Abtin Rahimian'
__license__   = '''
Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
'''

import argparse as ap
import logging
import os
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

lgr = logging.getLogger(__name__)

def parse_args():
    # one can use click, or cliff for better control
    p = ap.ArgumentParser(description='Loads and plots the shape dumped by ves3d code.')
    p.add_argument('-p', help='spherical harmonics order', type=int)
    p.add_argument('-n', help='number of surfaces', default=1, type=int)
    p.add_argument('--out-template', '-o', help='output file name', default=None)
    p.add_argument('--animate', '-a', help='animate the plot', action='store_true')
    p.add_argument('--increment', '-i', help='animate the plot', default=1, type=int)
    p.add_argument('files', nargs=ap.REMAINDER, help='list of files')

    args = p.parse_args()
    return vars(args)

def load_file(fnames):

    data = list()
    for f in fnames:
        data.append(np.genfromtxt(f, dtype=None).reshape((-1,1)))
    data = np.hstack(data)
    return data

def plot_series( data, p, n, out_template, animate,
                 increment, **kwargs):

    sz  = 2*p*(p+1)
    ds  = data.size
    nT  = int(ds/sz/n/4)
    nC  = data.shape[1]
    nR  = sz*n/nC

    fig = plt.figure()
    for iT in range(0,nT,increment):
        print('step %d/%d' % (iT, nT))
        ax=plt.gca(projection='3d')
        ax.cla()

        if animate:
            el=20-8*(iT/nT/2-1)**2
            az=120+160*iT/nT
            ax.view_init(elev=el,azim=az)

        d   = data[4*iT*nR:4*(iT+1)*nR,:]
        xyz = d[0:3*nR,:].reshape((3*sz,-1), order="F")
        s   = d[3*nR:,:].reshape((sz,-1), order="F")

        for iN in range(n):
            plot(p=p,
                 x = xyz[   0:  sz,iN],
                 y = xyz[  sz:2*sz,iN],
                 z = xyz[2*sz:    ,iN],
                 s = s[:,iN])

        plt.draw()
        ax.axis('tight')
        plt.pause(.05)

        if out_template:
            fname = out_template % iT
            plt.savefig(fname, transparent=True)

    plt.show()

def plot(p,x,y,z,s):

    x = x.reshape((p+1,-1))
    y = y.reshape((p+1,-1))
    z = z.reshape((p+1,-1))

    x = np.hstack((x, x[:,0,None]))
    y = np.hstack((y, y[:,0,None]))
    z = np.hstack((z, z[:,0,None]))

    ax = plt.gca(projection='3d')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, color='r',
                           linewidth=.1, antialiased=False)

def main():
    opts = parse_args()
    data = load_file(opts['files'])
    plot_series(data=data,**opts)

if __name__ == '__main__':
    main()
