#!/bin/bash
set -x
echo "Current date : $(date) @ $(hostname)"
echo "                                     "
echo "                                     "
echo "                                     "
echo "# arguments called with ---->  ${@}     "
echo "# \$1 ----------------------->  $1      "
echo "# \$2 ----------------------->  $2      "
echo "# path to me --------------->  ${0}     "
echo "# parent path -------------->  ${0%/*}  "
echo "# my name ------------------>  ${0##*/} "
echo "                                     "
echo "-------------------------------------"
echo "        read from userfile           "
echo "-------------------------------------"
localdir=$(pwd)
userfile=${localdir}/config.txt
#
#generate=$(sed '2q;d' $userfile)
testfile=$(sed '2q;d' $userfile)
login=$(sed '4q;d' $userfile)
hpcdir=$(sed '6q;d' $userfile)
qsubdir=$(sed '8q;d' $userfile)
localy2d=$(sed '10q;d' $userfile)
paraviewdir=$(sed '12q;d' $userfile)
giddir=$(which gid)
#gmshdir=$(which gmsh)
echo "                                     "
# get file name
file=$(basename -- "$testfile")
filename="${file%%.*}"
localresults="${localdir}/results/${filename}"
fullfilename="${localdir}/${filename}"
fulltestfile="${fullfilename}.y"
logfile="${fullfilename}.log"
#
rm $logfile                                  
echo "-------------------------------------"    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "            LOG FILE                 "    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "-------------------------------------"    >> $logfile
#rm -r $localresults                             >> $logfile
mkdir -p $localresults                          >> $logfile
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
    
#if [ $generate == 'yes' ]
#then
#    # execute gid: bad documentation, command unknown
#    #$localdir/gid/B2D.gid/GID_B2D.exe $fullfilename.dat $fullfilename.par >> $logfile
#    #
#    # generate the y input file
#    # g++ -o executable  -I/usr/include/gmsh
#    g++ $localdir/ygen.cpp $gmshdir
#fi
#
echo "                                     "    >> $logfile
# execute y2d locally?
if [ $localy2d == 'yes' ]
then
    # give permission to Y program
    chmod +x $localdir/Yf
    chmod +x $localdir/m2vtu
    chmod +x $localdir/m2vtu_crack
    echo "                      "
    echo "run job              "
    $localdir/Yf $fulltestfile
    echo "convert stress info into VTU files"
    $localdir/m2vtu $fullfilename $fulltestfile 1 100 1
    echo "convert crack info into VTU files"
    $localdir/m2vtu_crack $fullfilename $fulltestfile 1 100 1
echo "                                     "    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "            move files               "    >> $logfile
echo "-------------------------------------"    >> $logfile
echo "                                     "    >> $logfile
    #
    #
    #
    # move all files into the results folder
    #
    #
    #
    mv -v *.ym $localresults                            >> $logfile
    mv -v *.y $localresults                             >> $logfile
    mv -v .msh $localresults                            >> $logfile
    mv -v contactforce.txt $localresults                >> $logfile
    mv -v *.vtu $localresults                           >> $logfile
    mv -v *.vtu_crack $localresults                     >> $logfile
    mv -v Ytmp $localresults                            >> $logfile
else
    # replace every instance 'test' by '${filename}' in script file
    sed -i 's|$HOME/project|'${hpcdir}'|g' ${localdir}/script.sh    >> $logfile
    sed -i 's|localresults|'${localresults}'|g' ${localdir}/script.sh       >> $logfile
    #sed -i "s|test|${filename}|g" ${localdir}/scripts/script.sh                     >> $logfile
    #prepare hpc setup
    # create results folder
    hpcresults=${hpcdir}/results2d
    $login mkdir -p $hpcresults                             >> $logfile
        # get the setup from local directory
    $login rm -r $hpcdir                                    >> $logfile
    $login mkdir -p $hpcdir                                 >> $logfile
    scp -r $localdir $hpcdir                                >> $logfile
    ###########################################################################
    #run hpc
    ###########################################################################
    # need to call qsub with the path from hpc home
    shorthpcdir=${hpcdir##*/} #${hpcdir#*home}
    $login $qsubdir/qsub ${shorthpcdir}/script.sh          >> $logfile
    # revert changes to script file
    sed -i "s|${filename}|test|g" ${localdir}/script.sh                     >> $logfile
    sed -i 's|'${hpcdir}'|$HOME/project|g' ${localdir}/script.sh    >> $logfile
    sed -i 's|'${localresults}'|localresults|g' ${localdir}/script.sh       >> $logfile
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
