#!/bin/sh
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:mem=96gb

module load vtk/5.8.0
module load paramodule load paraview/5.4.1
mpiexec pvserver --reverse-connection --client-host=80.42.117

mkdir ${EPHEMERAL}/${PBS_JOBID}

cp -r $HOME/project/setup2d/* $EPHEMERAL/${PBS_JOBID}

userfile=user.txt
testfile=$(sed '2q;d' $userfile)
echo "$testfile"
login=$(sed '4q;d' $userfile)
echo "$login"
hpcdir=$(sed '6q;d' $userfile)
echo "$hpcdir"
localdir=$(sed '8q;d' $userfile)
echo "$localdir"
qsubdir=$(sed '10q;d' $userfile)
echo "$qsubdir"

testfilename=$(echo "$testfile" | cut -f 1 -d '.')
echo "$testfilename"

cd ${EPHEMERAL}/${PBS_JOBID}

chmod 777 Yf
chmod 777 m2vtu
chmod 777 m2vtu_crack

./Yf testfile
./m2vtu_crack testfilename testfile 1 100 1
./m2vtu testfilename testfile 1 100 1

mkdir -p $HOME/project/results2d/${PBS_JOBID}

cp ./* $HOME/project/results2d/${PBS_JOBID}

find $HOME/project/setup2d -type f -not -name 'script.sh' | xargs rm -rf