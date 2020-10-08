#!/bin/bash

gnuplot -e "
set log xy;
set key left bottom;
set xlabel 'Size [kB]';
set ylabel 'Operations per second';
set grid;
set terminal pngcairo;
set output 'results.png';
plot 'results.txt' i 0 u 2:(1e9*\$4) w l t 'Random permutation', 'results.txt' i 1 u 2:(1e9*\$4) w l t 'Sequential (4B)', 'results.txt' i 2 u 2:(1e9*\$4) w l t 'Sequential (64B)';
"
