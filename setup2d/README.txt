change variables names in the customisation file "custom.txt"

*******************************************************************************************
**************************** main script **************************************************
*******************************************************************************************
test.sh: bash script
- calculates y input file using Y2D code executables
	- either calculates in hpc directory, for this purpose, directories are copied to hpc and a job is submitted
	- or calculates locally, however, vtk 5.8 is necessary to be installed
- stores results
*******************************************************************************************
*******************************************************************************************

getresults.sh get results from the HPC system without having to log in manually

script.sh: bash script 
- submitted to the hpc system
- uses Y2D code to calculate the y file and create vtu files

qstat: bash script
- view status of submitted jobs without having to log into hpc manually
