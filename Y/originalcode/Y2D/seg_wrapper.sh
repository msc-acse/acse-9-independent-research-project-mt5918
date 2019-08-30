#!/bin/bash
gdbcommandfile=$(tempfile)
usepid=$(cat /proc/sys/kernel/core_uses_pid)
printf "set pagination off\nbacktrace\nquit\n" > $gdbcommandfile
ulimit -c unlimited
"$@"&
pid=$!
wait $!
if [[ $? -eq 139 ]]; then
    if [[ $usepid == 1 ]]; then 
        gdb -q $1 core.$pid -x $gdbcommandfile
    else
        gdb -q $1 core -x $gdbcommandfile
    fi
fi
rm $gdbcommandfile