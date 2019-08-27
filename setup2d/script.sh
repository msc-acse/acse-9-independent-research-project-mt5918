#!/bin/sh
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:mem=96gb
#PBS -N test

module load vtk/5.8.0
module load paramodule load paraview/5.6.1

mkdir ${EPHEMERAL}/${PBS_JOBID}

cp $HOME/project/setup2d/* $EPHEMERAL/${PBS_JOBID}

cd ${EPHEMERAL}/${PBS_JOBID}

chmod 777 Yf
chmod 777 m2vtu
chmod 777 m2vtu_crack

./Yf test.y
./m2vtu test test.y 1 100 1
./m2vtu_crack test test.y 1 100 1

mkdir -p $HOME/project/setup2d/../results

cp ./* $HOME/project/setup2d/../results

#while [ 1 ]; do
#CLIENT=$(cat $HOME/pvclient.txt)
#pvserver --force-offscreen-rendering -rc --client-host=$CLIENT --sp=1111
#done
