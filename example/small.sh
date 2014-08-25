#!/bin/bash -x
#
#PJM --rsc-list "node=24"
#PJM --rsc-list "elapse=00:10:00"
#PJM --rsc-list "rscgrp=small"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin "rank=* ../smash %r:./"
#PJM --stgin "rank=* ./small.inp %r:./"
#PJM -s
#PJM --mpi "proc=24"

. /work/system/Env_base

export OMP_NUM_THREADS=8

mpiexec -stdin small.inp lpgparm -s 4MB -d 4MB -h 4MB -t 4MB -p 4MB ./smash

