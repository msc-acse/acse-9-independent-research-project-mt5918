#!/bin/bash
userfile=config.txt
login=$(sed '4q;d' $userfile)
qsubdir=$(sed '8q;d' $userfile)
$login $qsubdir/qstat
sleep 10s
