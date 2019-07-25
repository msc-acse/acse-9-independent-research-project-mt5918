#!/bin/bash
userfile=../user/config.txt
login=$(sed '6q;d' $userfile)
qsubdir=$(sed '12q;d' $userfile)
$login $qsubdir/qstat
sleep 10s
