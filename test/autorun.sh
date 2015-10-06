#!/bin/bash
#
# check the designated file and run it when it changes.

if [ "$#" -eq 1 ]; then
    check=$1
    run=$1
else
    check=$1
    run="${*:2}"
fi

# try -m kqueue_monitor for fswatch if the default is not working
fswatch -to -l 0.2 -E ${check} | while read e; do clear; echo "running '${run}' on" `date`; ${run}; done
