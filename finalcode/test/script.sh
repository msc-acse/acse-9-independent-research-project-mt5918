#!/bin/sh

cd ../Y2D
make clean
make

read

mkdir -p ../test/results
cd ../test/results 

cp ../../Y2D/* .

chmod +x Yf

#clear

./Yf test.y

find . -type f ! -name "test*" -delete