#!/bin/bash
#$ -t 1-100
#$ -N lorenz96_cpp
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -m beas
#$ -o ./output
#$ -e ./error
#$ -q batch.q
./lorenz96_cpp $SGE_TASK_ID
