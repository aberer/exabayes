#! /bin/sh 

alts="mpirun mpiexec mpiexec-openmpi-mp  mpiexec-mpich-mp"

for elem in `echo $alts `
do 
    loc=`which $elem`
    if [ "$loc" != "" ]; then
	echo $loc
	exit 
    fi
done 

echo "Could not find mpirun or mpiexec. Cannot run parallel programs, unless mpi is installed." >&2
exit 1 

