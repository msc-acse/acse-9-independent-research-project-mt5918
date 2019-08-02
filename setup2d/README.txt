this directory is: {localdirectory}/scripts

variables names in curly brackets { } can be changed in the customisation file "custom.txt"

*******************************************************************************************
**************************** main script **************************************************
*******************************************************************************************
test.sh: bash script
- (optional) generates a y testfile in {localdirectory}/gid/{testfile}/{testfile}.y 
- calculates {localdirectory}/gid/{testfile}/{testfile}.y using Y2D code executables from {localdirectory}/bin
	- either calculates in hpc directory, for this purpose, directories are copied to hpc and a job is submitted
	- or calculates locally, however, vtk 5.8 is necessary to be installed
- stores results in {localdirectory}/gid/{testfile}/results
- visualises results with paraview
*******************************************************************************************
*******************************************************************************************

test.bat: corresponding windows script
- not working properly // abandonned

script.sh: bash script 
- submitted to the hpc system
- uses Y2D code to calculate the y file and create vtu files

qsub: bash script
- submit a local job directly to the hpc system without having to log into hpc manually

qstat: bash script
- view status of submitted jobs without having to log into hpc manually
