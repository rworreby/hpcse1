#!/usr/bin/env python
# File       : measure_speedup.py
# Description: Create a pdf plot from measurements collected with the
#              measure_speedup.sh BASH script.
# Copyright 2020 ETH Zurich. All Rights Reserved.
import argparse
import glob
import numpy as np
import matplotlib as mpl
mpl.use('Agg') # render without X server running
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-m', '--measurements', required=True, nargs='+',
            help='Individual measurement files')
    parser.add_argument('-p', '--prefix', type=str, default='',
            help='Prefix used for output files')

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    test_size = ['small', 'large']
    lanes     = ['2-way', '4-way']
    results   = {
            'threads': [],
            'small': { # first test size
                '2-way':{'mean': [], 'std':[]},
                '4-way':{'mean': [], 'std':[]}
                },
            'large': { # second test size
                '2-way':{'mean': [], 'std':[]},
                '4-way':{'mean': [], 'std':[]}
                }
            }
    for measurement in sorted(args.measurements):
        results['threads'].append( int(measurement.split('_')[1]) )
        data = np.loadtxt(measurement)
        for i,size in enumerate(test_size):
            for j,simd in enumerate(lanes):
                results[size][simd]['mean'].append( np.mean(data[:,2*i+j]) )
                results[size][simd]['std'].append( np.std(data[:,2*i+j]) )

    # plot
    cm = plt.cm.coolwarm
    cidx = np.linspace(0.0, 1.0, len(lanes))
    fmt = ['-o', '-s']
    for simd in lanes:
        SIMD = float(simd.split('-')[0])
        tmin = float(min(results['threads']))
        tmax = float(max(results['threads']))
        for i,size in enumerate(test_size):
            plt.errorbar(results['threads'], results[size][simd]['mean'],
                    yerr=results[size][simd]['std'], fmt=fmt[i], color=cm(cidx[i]),
                    lw=1.5, ms=6, capsize=3, capthick=1.5,
                    label='Test array size {:s}'.format(size))
        plt.plot([tmin, tmax], SIMD*np.array([tmin, tmax]), color='k', lw=1.5, ls='--', label='Ideal')
        plt.grid(which='both')
        plt.xlabel('Number of Threads')
        plt.ylabel('Speedup')
        plt.legend(loc='upper left')
        plt.savefig('{:s}speedup_{:s}.pdf'.format(args.prefix, simd), bbox_inches='tight')
        plt.close()
