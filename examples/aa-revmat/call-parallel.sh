#! /bin/bash

exec=""

if [ -d ../../build -a -f ../../build/exabayes ]; then    
    exec="../../build/exabayes"
elif [ -f ../../exabayes ]; then
    exec="../../exabayes"
else 
    topname=$(cd ../../ ; pwd )
    echo <<EOF "Error: could not find exabayes in either $topname or
$topname/build.  Please use the build script to build the sequential
version of ExaBayes or build it in the top-level exabayes-directory."
EOF

    exit
fi


mpi=""
if [ "$(which mpirun)"  != "" ]; then
    mpi=mpirun
elif [ "$(which openmpirun)" != "" ] ; then 
    mpi=openmpirun
else 
    echo <<EOF
Error: Could not find mpi. 
EOF
    exit 
fi


# this is a pretty memory intensive dataset. Let's trade as much runtime for memory as possible: 

cmd="$mpi -np 2 $exec -f aln.phy -m PROT -c config.nex -n myRun -s 123 -M 3 -S"
eval $cmd
