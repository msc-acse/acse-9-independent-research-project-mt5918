#!/bin/sh

set -e

cd ../Y2D
make clean
make

read

mkdir -p ../test/results
cd ../test/results 

cp ../../Y2D/* .

./Yf test.y