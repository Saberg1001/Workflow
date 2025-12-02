#!/bin/bash
#SBATCH -J gaussian-task 
#SBATCH -p XDD           
#SBATCH -N 1            
#SBATCH -c 8             

export GAUSS_EXEDIR=/data/gauss/g16 
export GAUSS_SCRDIR=$PWD 

g16 jobname output.log  