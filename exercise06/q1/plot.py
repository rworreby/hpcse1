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
    u = np.fromfile(path, dtype=np.float64)
    N = int(np.round((len(u) + 1)**0.5 - 1))
    fullN = N + 2
    assert len(u) == N * fullN, \
            "Couldn't deduce shape: got {:} elements, expected {:} for N={:}".\
            format(len(u), N * fullN, N)
    u = u.reshape((N, fullN))
    plt.imshow(np.flipud(u.T), extent=(0, 1, 0, fullN / N), vmin=0, vmax=1)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.axis('off')
    plt.axis('equal')
    fpath = stripped_path + '.png'
    print(fpath)
    plt.savefig(fpath, bbox_inches='tight')
    plt.close()


def main():
    args, _ = parseArgs()

    for f in args.files:
        plot_dat(f)


if __name__ == '__main__':
    main()
