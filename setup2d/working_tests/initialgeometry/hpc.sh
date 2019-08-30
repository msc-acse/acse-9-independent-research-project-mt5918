#!/bin/sh

module load vtk/5.8.0
module load paramodule load paraview/5.6.1

chmod 777 Yf            
chmod 777 m2vtu         
chmod 777 m2vtu_crack   

./Yf test3                           
./m2vtu test3 test3.y 1 100 1         
./m2vtu_crack test3 test3.y 1 100 1   

mkdir -p results2d/test3
cp * results2d/

#while [ 1 ]; do
#CLIENT=$(cat $HOME/pvclient.txt)
#pvserver --force-offscreen-rendering -rc --client-host=$CLIENT --sp=1111
#done
