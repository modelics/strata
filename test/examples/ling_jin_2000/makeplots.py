#!/usr/bin/env python

# Author: Shashwat Sharma

# Copyright 2021 Shashwat Sharma and Piero Triverio

# This file is part of Strata.

# Strata is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Strata is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Strata.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import csv
import os
from math import pi


def PlotAllComponents(refname, reffiles, testfile, refscales, component_idx, component_names):

    rc("text", usetex = True)

    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 16}

    rc('font', **font)
    
    eps0 = 8.854*1e-12
    mu0 = 4.0*pi*1e-7

    with open(testfile) as temp:
        f = temp.readline()
        f = float(f.split(" ")[1:][0])

    k0 = 2.0*pi*f*np.sqrt(eps0*mu0)

    fig, ax = plt.subplots(figsize=(8, 4))

    idx = -1
    for reffile in reffiles:

        idx = idx + 1

        name, extension = os.path.splitext(reffile)

        delim = ' '
        if extension == ".csv":
            delim = ','

        refdata = np.loadtxt(reffile, delimiter=',')
        testdata = np.loadtxt(testfile, delimiter=' ', skiprows=5)

        x_ref = np.asarray(refdata[:,0])
        y_ref = np.asarray(refdata[:,1])
        y_ref = y_ref*refscales[idx]

        x_test = np.asarray(testdata[:,0])
        x_test = x_test*k0
        y_test = np.asarray(testdata[:,component_idx[idx]])
        
        if idx == 0:
            reflegend = refname
            testlegend = 'Strata'
        else:
            reflegend = ''
            testlegend = ''

        ax.plot(x_ref, y_ref, 'k-', label=reflegend)
        ax.plot(x_test, y_test, 'r--', label=testlegend)

    ax.set(xlabel = "$k_0\\rho$",
           ylabel = "$|G|$",
           yscale = 'log',
           xscale = 'log'
    )

    ax.grid(b=True, which='both', axis='both', linewidth=0.4)
    plt.legend(loc='best', prop={'size': 14})
    
    plt.tight_layout()
    fig.savefig("./plots.pdf", bbox_inches='tight')
    
    return


def main():
    
    refname = "Ling, Jin, 2000, Fig. 3"

    reffiles = ["./fig3_Gxx.csv", "./fig3_Gzx.csv", "./fig3_Gzz.csv", "./fig3_Gphi.csv"]
    testfile = "./MGFdata.txt"
    # testfile = "./MGFdata_interp.txt"
    # testfile = "./MGFdata_DCIM.txt"
    # testfile = "./MGFdata_singularity.txt"

    refscales = [1.0, 1.0, 1.0, 1.0]

    component_idx = [1, 7, 9, 10];
    component_names = ["G_{xx}", "G_{zx}", "G_{zz}", "G_{\\phi}"]

    PlotAllComponents(refname, reffiles, testfile, refscales, component_idx, component_names)

    plt.show()

    return


if __name__ == "__main__":
    main()

