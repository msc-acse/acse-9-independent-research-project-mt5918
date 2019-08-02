#!/bin/bash
set -x

# get file name
filename=$(basename -- "$fullfile")
filename="${filename%%.*}"
logfile="$filename.log"

#create log file
rm logfile
rm ../1.tmp

echo "-------------------------------------" >> $logfile
echo "        read from userfile           " >> $logfile
echo "-------------------------------------" >> $logfile
userfile=../../user/config.txt >> $logfile

generate=$(sed '2q;d' $userfile)
testfile=$(sed '4q;d' $userfile)
login=$(sed '6q;d' $userfile)
hpcdir=$(sed '8q;d' $userfile)
localdir=$(sed '10q;d' $userfile)
qsubdir=$(sed '12q;d' $userfile)
bindir=$(sed '16q;d' $userfile)
mayavidir=$(sed '20q;d' $userfile)

echo "testfile: $testfile" >> $logfile
echo "login: $login" >> $logfile
echo "hpcdir: $hpcdir" >> $logfile
echo "localdir: $localdir" >> $logfile
echo "qsubdir: $qsubdir" >> $logfile
echo "bindir: $bindir" >> $logfile
echo "mayavidir: $mayavidir" >> $logfile

echo "-------------------------------------" >> $logfile
echo "            run program              " >> $logfile
echo "-------------------------------------" >> $logfile

$filename.exe $filename.dat $filename.par

if exist 1.tmp goto end

chmod +x ./../../bin/Y
chmod +x ./../../bin/m2vtu2D
chmod +x ./../../bin/m2vtu2D_crack

echo "run job              " >> $logfile
. ./../../bin/Y $testfile
echo "convert stress info into VTU files" >> $logfile
. ./../../bin/m2vtu2D $filename $testfile
echo "convert crack info into VTU files" >> $logfile
. ./../../bin/m2vtu2D_crack $filename $testfile
echo "visualize results" >> $logfile
"$mayavidir/mayavi2" -d ${filename}0.vtu -m SurfaceMap -f ExtractTensorComponents -d ${filename}_crack0.vtu -m SurfaceMap
:end

