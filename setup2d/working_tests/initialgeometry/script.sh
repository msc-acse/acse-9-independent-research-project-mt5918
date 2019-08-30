#!/bin/sh
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:mem=96gb
#PBS -N test3
###PBS -J 1-10
### $PBS_ARRAY_INDEX

module load vtk/5.8.0
module load paramodule load paraview/5.6.1

mkdir -p ${EPHEMERAL}/test3 

cp -r ./* $EPHEMERAL/test3

cd ${EPHEMERAL}/test3 

chmod 777 Yf            
chmod 777 m2vtu         
chmod 777 m2vtu_crack   

./Yf test3                           
./m2vtu test3 test3.y 1 100 1         
./m2vtu_crack test3 test3.y 1 100 1   

mkdir -p $HOME/hpc/results2d

cp ./* $HOME/hpc/results2d

#while [ 1 ]; do
#CLIENT=$(cat $HOME/pvclient.txt)
#pvserver --force-offscreen-rendering -rc --client-host=$CLIENT --sp=1111
#done
