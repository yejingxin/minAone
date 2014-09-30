#!/bin/bash
#$ -t 1-100
#$ -N lorenz96-m1
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -m beas
#$ -o ./output
#$ -e ./error
#$ -q batch.q
./lorenz96_cpp $SGE_TASK_ID
