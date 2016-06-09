#!/bin/bash
#$ -t 1-20
#$ -N hh-m1
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -m beas
#$ -o ./output
#$ -e ./error
#$ -q batch.q
./simple_nakl_cpp $SGE_TASK_ID
