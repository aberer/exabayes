#! /bin/bash

if [ "$#" != "1" ]; then
    echo "$0 ompi|mpich"
    exit 
fi


mode=$1

case $mode in
    ompi)
	sudo update-alternatives --set mpirun /usr/bin/mpirun.openmpi  
	sudo update-alternatives --set mpi /usr/lib/openmpi/include
	;;
    mpich)
	sudo update-alternatives --set mpirun /usr/bin/mpirun.mpich
	sudo update-alternatives --set mpi /usr/include/mpich
	;;
    *)
	echo "$0 opmpi|mpich"
	exit
esac
