#!/bin/bash
set -x
echo "Current date : $(date) @ $(hostname)"
echo "                                     "
echo "                                     "
echo "                                     "
echo "# arguments called with ---->  ${@}     "
echo "# \$1 ---------------------->  $1       "
echo "# \$2 ---------------------->  $2       "
echo "# path to me --------------->  ${0}     "
echo "# parent path -------------->  ${0%/*}  "
echo "# my name ------------------>  ${0##*/} "
echo "                                     "
echo "-------------------------------------"
echo "        read from userfile           "
echo "-------------------------------------"
cd ..
localdir=$(pwd)
userfile=${localdir}/user/config.txt
#
generate=$(sed '2q;d' $userfile)
testfile=$(sed '4q;d' $userfile)
login=$(sed '6q;d' $userfile)
hpcdir=$(sed '8q;d' $userfile)
qsubdir=$(sed '12q;d' $userfile)
mayavidir=$(sed '16q;d' $userfile)
localy2d=$(sed '18q;d' $userfile)
giddir=$(sed '22q;d' $userfile)
y2ddir=$(sed '26q;d' $userfile)
paraviewdir=$(sed '30q;d' $userfile)
echo "                                     "
# get file name
file=$(basename -- "$testfile")
filename="${file%%.*}"
localresults="${localdir}/gid/${filename}.gid/results"
fullfilename="${localdir}/gid/${filename}.gid/${filename}"
fulltestfile="${fullfilename}.y"   
logfile="${fullfilename}.unix.log"
rm $logfile                                  
echo "-------------------------------------"    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "            LOG FILE                 "    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "-------------------------------------"    >> $logfile
rm -r $localresults                             >> $logfile
mkdir -p $localresults                             >> $logfile
cd $localresults                                >> $logfile
echo "                                     "    >> $logfile
echo "                                     "    >> $logfile
echo "                                     "    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "            Variables                "    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "testfile: $testfile"                      >> $logfile
echo "login: $login"                            >> $logfile
echo "hpcdir: $hpcdir"                          >> $logfile
echo "localdir: $localdir"                      >> $logfile
echo "qsubdir: $qsubdir"                        >> $logfile
echo "mayavidir: $mayavidir"                    >> $logfile
echo "                                     "    >> $logfile
echo "file: $file"                              >> $logfile
echo "filename: $filename"                      >> $logfile
echo "localresults: $localresults"              >> $logfile
echo "fullfilename: $fullfilename"              >> $logfile
echo "fulltestfile: $fulltestfile"              >> $logfile
echo "logfile: $logfile"                        >> $logfile
echo "                                     "    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "            run program              "    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "                                     "    >> $logfile
    
if [ $generate == 'yes' ]
then
    $localdir/gid/B2D.gid/GID_B2D.exe $fullfilename.dat $fullfilename.par >> $logfile
fi
echo "                                     "    >> $logfile
if [ $localy2d == 'yes' ]
then
    chmod +x $localdir/bin/Yf
    chmod +x $localdir/bin/m2vtu
    chmod +x $localdir/bin/m2vtu_crack
    echo "               s      "
    echo "run job              "
    $localdir/bin/Yf $fulltestfile
    echo "convert stress info into VTU files"
    $localdir/bin/m2vtu $fullfilename $fulltestfile 1 100 1
    echo "convert crack info into VTU files"
    $localdir/bin/m2vtu_crack $fullfilename $fulltestfile 1 100 1
echo "                                     "    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "            move files               "    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "                                     "    >> $logfile
    mv -v ../*.ym $localresults                            >> $logfile
    mv -v ../*.y $localresults                             >> $logfile
    mv -v ../.msh $localresults                            >> $logfile
    mv -v ../contactforce.txt $localresults                >> $logfile
    mv -v ../*.vtu $localresults                           >> $logfile
    mv -v ../*.vtu_crack $localresults                     >> $logfile
    mv -v ../Ytmp $localresults                            >> $logfile
else
    #replace every instance 'test' by '${filename}' in script file
    sed -i "s|test|${filename}|g" ${localdir}/scripts/script.sh                     >> $logfile
    sed -i 's|$HOME/project|'${hpcdir}'|g' ${localdir}/scripts/script.sh    >> $logfile
    sed -i 's|localresults|'${localresults}'|g' ${localdir}/scripts/script.sh       >> $logfile
    hpcresults=${hpcdir}/results2d
    $login mkdir -p $hpcresults                             >> $logfile
    #prepare hpc setup
    $login rm -r $hpcdir/setup2d                            >> $logfile
    $login mkdir -p $hpcdir/setup2d                                 >> $logfile
    scp -r $localdir $hpcdir/setup2d                                >> $logfile
    #run hpc
    shorthpcdir=${hpcdir#*home}
    $login $qsubdir/qsub project/setup2d/scripts/script.sh          >> $logfile
    
    sed -i "s|${filename}|test|g" ${localdir}/scripts/script.sh                     >> $logfile
    sed -i 's|'${hpcdir}'|$HOME/project|g' ${localdir}/scripts/script.sh    >> $logfile
    sed -i 's|'${localresults}'|localresults|g' ${localdir}/scripts/script.sh       >> $logfile
fi
echo "                                     " >> $logfile
echo "-------------------------------------" >> $logfile
echo "            visualize results        " >> $logfile
echo "-------------------------------------" >> $logfile
echo "                                     " >> $logfile
#wait
#"$mayavidir/mayavi2" -d ${fullfilename}0.vtu -m SurfaceMap -f ExtractTensorComponents -d ${fullfilename}_crack0.vtu -m SurfaceMap   >> $logfile
#"$paraviewdir/paraview.exe" $localresults/${filename}*.vtu $localresults/${filename}_crack*.vtu      >> $logfile
echo "                                     "            >> $logfile
echo "Done                                 "            >> $logfile
