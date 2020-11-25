#!/usr/bin/env bash
# File       : measure_speedup.sh
# Description: Measure vectorized speedup with multiple cores.  Once you have
#              completed your code, you can run this script with
#              $ make measurement
#              using the 'measurement' target in the Makefile.  It will
#              schedule a job on a compute node on euler.  See the Makefile.
# Copyright 2020 ETH Zurich. All Rights Reserved.
export OMP_DYNAMIC='false'
export OMP_PROC_BIND='true'

get_sample()
{
    result="${1}"; shift
    sample="/tmp/sample${RANDOM}.out"

    sleep 3 # creates some space between measurements
    ./vec_red | tee ${sample} # run vectorized reduction code

    small_2way="$(awk 'NR==26 {print $2}' ${sample})"
    small_4way="$(awk 'NR==20 {print $2}' ${sample})"
    large_2way="$(awk 'NR==53 {print $2}' ${sample})"
    large_4way="$(awk 'NR==47 {print $2}' ${sample})"

    rm -f ${sample}
    echo "${small_2way} ${small_4way} ${large_2way} ${large_4way}" >> ${result}
}

for threads in 1 2 4 8 12 16 20 24
do
    export OMP_NUM_THREADS=${threads}
    printf -v result "/tmp/threads${LSB_JOBID}_%03d" "${threads}"
    for sample in $(seq 10) # collect 10 samples (~30 seconds to collect all)
    do
        get_sample ${result}
    done
done

# make plot
module load new
module load python/3.6.1
python plot_speedup.py \
    --prefix "${LSB_JOBID}" \
    --measurements /tmp/threads${LSB_JOBID}_* && \
    rm -f /tmp/threads${LSB_JOBID}_* # clean up
