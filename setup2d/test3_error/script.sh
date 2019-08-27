#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=96gb
#PBS -N test3
module load vtk/5.8.0
mkdir ${EPHEMERAL}/${PBS_JOBNAME}
cp -r $HOME/test3/* $EPHEMERAL/${PBS_JOBNAME}
cd ${EPHEMERAL}/${PBS_JOBNAME}
chmod 777 Yf
chmod 777 m2vtu
chmod 777 m2vtu_crack
./Yf test3.y
./m2vtu_crack test3 test3.y 1 100 1
./m2vtu test3 test3.y 1 100 1
mkdir -p $HOME/FEMDEM_2d_REs/${PBS_JOBNAME}
cp ./* $HOME/FEMDEM_2d_REs/${PBS_JOBNAME}
