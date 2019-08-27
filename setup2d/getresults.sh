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
localdir=$(pwd)
userfile=${localdir}/config.txt
#
#generate=$(sed '2q;d' $userfile)
testfile=$(sed '2q;d' $userfile)
login=$(sed '4q;d' $userfile)
hpcdir=$(sed '6q;d' $userfile)
qsubdir=$(sed '8q;d' $userfile)
giddir=$(which gid)
paraviewdir=$(sed '12q;d' $userfile)
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
scp -r $hpcdir/results $localdir                >> $logfile