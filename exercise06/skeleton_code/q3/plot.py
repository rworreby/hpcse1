#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
import os


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('files',
                        type=str,
                        nargs='+',
                        help='Paths to files dump_*.dat')
    return parser.parse_known_args()


def plot_dat(path):
    stripped_path = os.path.splitext(path)[0]
    with open(path, 'r') as header:
        t = float(header.readline().strip().split()[2])
        s = int(header.readline().strip().split()[2])
        N = int(header.readline().strip().split()[2])
    f = np.loadtxt(path).reshape((N, N))
    cs = plt.imshow(np.flipud(f), extent=(0, 1, 0, 1), vmin=0, vmax=1)
    plt.colorbar(cs)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.axis('off')
    plt.axis('equal')
    plt.title(f'Time: {t}')
    fpath = stripped_path + '.png'
    print(fpath)
    plt.savefig(fpath, dpi=200, bbox_inches='tight')
    plt.close()


def main():
    args, _ = parseArgs()

    for f in args.files:
        plot_dat(f)


if __name__ == '__main__':
    main()
