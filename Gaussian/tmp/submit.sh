#!/bin/bash
#SBATCH -J gaussian-task
#SBATCH -e err
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p XDM004

export g16root=/data/public-software/gauss/
source $g16root/g16/bsd/g16.profile

echo "USER: $USER"
echo "SLURM_JOBID: $SLURM_JOBID"

export GAUSS_SCRDIR=/tmp/$USER/$SLURM_JOBID
mkdir -p "$GAUSS_SCRDIR"

INPUTDIR=`pwd`
USERNAME=`whoami`
DATA_TAG=`date +%F_%T|sed 's/:/./g'`


# run g16
g16 jobname output.log

formchk neutral.chk || true
formchk anion.chk   || true
formchk cation.chk  || true
formchk excited.chk || true

nt_count=$(grep -c "Normal termination of Gaussian" output.log || echo 0)
has_conv_fail=$(grep -q "Convergence criterion not met" output.log && echo 1 || echo 0)

record_line="$DATA_TAG    Gaussian job ${SLURM_JOBID:-N/A}    $INPUTDIR"

if [ "$nt_count" -eq 5 ] && [ "$has_conv_fail" -eq 0 ]; then
    echo "$record_line" >> ~/Gaussian/finish.log
else
    echo "$record_line" >> ~/Gaussian/unfinish.log
fi

