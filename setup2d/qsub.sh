#!/bin/bash
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
scp -r $localdir $hpcdir
$login $qsubdir/qsub project/setup2d/script.sh