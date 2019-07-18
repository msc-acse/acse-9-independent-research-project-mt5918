#!/bin/bash
userfile=user.txt
login=$(sed '4q;d' $userfile)
qsubdir=$(sed '10q;d' $userfile)
$login $qsubdir/qstat
sleep 10s