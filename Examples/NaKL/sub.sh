#!/bin/bash
#$ -t 1-100
#$ -N hh-m1
#$ -cwd
#$ -j y
#$ -M yejingxin.ucsd@gmail.com
#$ -S /bin/bash
#$ -m beas
#$ -o ./output
#$ -e ./error
#$ -q batch.q
./simple_nakl_cpp $SGE_TASK_ID
