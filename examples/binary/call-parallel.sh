#! /bin/bash

# find the executable 
dir=`dirname $0`
bin=`$dir/../findBin.sh para`
if [ $? != 0  ]; then
    exit
fi 

# find mpi run 
mpi=`$dir/../findMpi.sh`
if [ $? != 0 ]; then
    exit 
fi


# this is a pretty memory intensive dataset. Let's trade as much runtime for memory as possible: 
cmd="$mpi -np 8  $bin -R 2 -C 2  -f aln.phy -q aln.part -c config.nex -n myRun -s 123"
echo  "\n\ncommandline:\n$cmd\n\n"
eval $cmd
