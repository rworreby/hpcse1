#!/bin/sh

set -eu

N=256
T=0.05

for dt in 0.01 0.005 0.0025 0.00125 ; do
  echo dt=$dt
  c="./diffusionADI  $N  $dt  $T  0"
  echo "$c"
  eval "$c"

  mv stat.dat stat_${dt}.dat
  echo
done

python3 plot_stat.py stat_*.dat
