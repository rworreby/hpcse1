#!/bin/sh

set -eu

function check()
{
  echo "Verifying results..."
  f=dump_sequential.dat
  g=dump_parallel.dat
  if ! diff "$f" "$g" ; then
    exit 1
  fi
  echo "Verification Passed"
}

for N in $(seq 1024 512 6144) ; do
  cmd="mpirun -n 8 ./diffusion 1.0 2.0 $N ; check"
  echo "$cmd"
  eval "$cmd"
  echo
done
