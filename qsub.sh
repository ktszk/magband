#!/bin/sh
#$ -N suscep3d
#$ -N job_name
#$ -pe smp 40
#$ -cwd
#$ -V
#$ -q all.q
#$ -S /bin/bash

export OMP_NUM_THREADS=$NSLOTS
#export OMP_NUM_THREADS=32
./main