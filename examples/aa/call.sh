#! /bin/bash

exec=""

if [ -d ../../build -a -f ../../build/yggdrasil ]; then    
    exec="../../build/yggdrasil"
elif [ -f ../../yggdrasil ]; then
    exec="../../yggdrasil"
else 
    topname=$(cd ../../ ; pwd )
    echo <<EOF "Error: could not find yggdrasil in either $topname or
$topname/build.  Please use the build script to build the sequential
version of ExaBayes or build it in the top-level exabayes-directory."
EOF

    exit
fi

# this is a pretty memory intensive dataset. Let's trade as much runtime for memory as possible: 

cmd="$exec -f aln.phy -m PROT -c config.nex -n myRun -s 123 -M 3 -S"
eval $cmd
