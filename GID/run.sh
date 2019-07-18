#!/bin/bash
shopt -s expand_aliases
alias icl='ssh mt5918@login.cx1.hpc.ic.ac.uk'
icl rm project/setup2d/test.y
scp /GID/test.y mt5918@login.cx1.hpc.ic.ac.uk:/rds/general/user/mt5918/home/project/setup2d
icl /opt/pbs/bin/qsub /rds/general/user/mt5918/home/project/setup2d/script.sh
wait
