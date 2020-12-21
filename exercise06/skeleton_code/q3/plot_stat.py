#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('files',
                        type=str,
                        nargs='+',
                        help='Paths to files stat_*.dat')
    parser.add_argument('--output',
                        type=str,
                        default="stat.png",
                        help='Output png path')
    return parser.parse_known_args()


def main():
    args, _ = parseArgs()

    for f in args.files:
        u = np.genfromtxt(f, names=True)
        t = u['t']
        dt = t[1] - t[0]
        plt.plot(t, u['max'], label='dt = ' + str(dt))

    plt.legend()
    plt.ylim([0, None])

    print(args.output)
    plt.savefig(args.output, dpi=200, bbox_inches='tight')


if __name__ == '__main__':
    main()
