#!/bin/sh
module load vtk/5.8.0
mkdir $HOME/FEMDEM_2d_REs/test
rm $HOME/FEMDEM_2d_REs/test/*
cp ./* $HOME/FEMDEM_2d_REs/test
cd $HOME/FEMDEM_2d_REs/test
chmod +x Yf
chmod +x m2vtu
chmod +x m2vtu_crack
./Yf test3.y
./m2vtu_crack test3 test3.y 1 100 1
./m2vtu test3 test3.y 1 100 1
rm $HOME/FEMDEM_2d_REs/{Yf,m2vtu,m2vtu_crack,*.y}
