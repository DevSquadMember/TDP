#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage : $0 <nom du fichier>"
else
  rm -f "err/$1.err"
  rm -f "out/$1.out"
  mpicc "$1.c" -o "bin/$1"
  sbatch "$1.slurm"
fi

