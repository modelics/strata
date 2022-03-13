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

    # In order to make the plots correspond to Fig. 3 and Fig. 5 in Aksun 1996
    tlim = [10, 20]

    idx = -1
    for reffile in reffiles:

        idx = idx + 1

        name, extension = os.path.splitext(reffile)

        delim = ' '
        if extension == ".csv":
            delim = ','

        # refdata = np.loadtxt(reffile, delimiter=',')
        testdata = np.loadtxt(testfile, delimiter=' ', skiprows=5)

        # x_ref = np.asarray(refdata[:,0])
        # y_ref = np.asarray(refdata[:,1])
        # y_ref = y_ref*refscales[idx]

        x_test = np.asarray(testdata[:,0])
        y_real_test = np.asarray(testdata[:,component_idx[idx]])
        y_imag_test = np.asarray(testdata[:,component_idx[idx]+1])

        fig, ax = plt.subplots(figsize=(8, 4))
   
        reflegend = refname
        testlegend = 'Strata'

        # ax.plot(x_ref, y_ref, 'k-', label=reflegend)
        ax.plot(x_test, y_real_test, 'b--', label=testlegend + ", real part")
        ax.plot(x_test, y_imag_test, 'r--', label=testlegend + ", imag. part")

        ax.set(xlabel = "$t$",
               ylabel = "$G$",
               yscale = 'linear',
               xscale = 'linear'
        )

        ax.grid(b=True, which='both', axis='both', linewidth=0.4)
        plt.legend(loc='best', prop={'size': 14})
    
        plt.tight_layout()
        fig.savefig("./plots_spectral.pdf", bbox_inches='tight')
    
    return


def main():
    
    refname = "Yuan, Sarkar, Salazar-Palma, 2006, Fig. 3"

    # reffiles = ["./fig3_Gxx.csv", "./fig3_Gzx.csv", "./fig3_Gzz.csv", "./fig3_Gphi.csv"]
    reffiles = ["temp"]
    testfile = "./MGFdata_spectral.txt"

    refscales = [1.0, 1.0, 1.0, 1.0]

    component_idx = [1];
    component_names = ["Gphi_Fig3"]

    PlotAllComponents(refname, reffiles, testfile, refscales, component_idx, component_names)

    plt.show()

    return


if __name__ == "__main__":
    main()

