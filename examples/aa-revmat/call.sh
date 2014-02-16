#! /bin/sh

# find the executable 
dir=`dirname $0`
bin=`$dir/../findBin.sh seq`
if [ $? != 0  ]; then
    exit
fi

# this is a pretty memory intensive dataset. Let's trade as much runtime for memory as possible: 
cmd="$bin -f aln.phy -m PROT -c config.nex -n myRun -s 123 -M 3 -S"
echo  "\n\ncommandline:\n$cmd\n\n"
eval $cmd
