#!/bin/sh
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:mem=96gb
#PBS -N test

module load vtk/5.8.0
module load paramodule load paraview/5.6.1

mkdir -p ${EPHEMERAL}/test 

cp -r $HOME/project/setup2d/bin $EPHEMERAL/test
cp -r $HOME/project/setup2d/gid/test.gid $EPHEMERAL/test

cd ${EPHEMERAL}/test 

chmod 777 Yf            
chmod 777 m2vtu         
chmod 777 m2vtu_crack   

./Yf test                           
./m2vtu test test.y 1 100 1         
./m2vtu_crack test test.y 1 100 1   

mkdir -p $HOME/project/results2d/test

cp ./* $HOME/project/results2d/test  

find $HOME/project/setup2d rm -r *   

#while [ 1 ]; do
#CLIENT=$(cat $HOME/pvclient.txt)
#pvserver --force-offscreen-rendering -rc --client-host=$CLIENT --sp=1111
#done
