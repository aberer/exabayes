#! /bin/bash


if [ $# -lt 3 ]; then
    echo -e  "$0 isApple cores ssetypes..\n\nwhere ssetypes can be  avx,sse, no-sse"
    exit
fi

IS_APPLE=$1
shift
cores=$1
shift

ssetypes=$*

fold=$(pwd  | tr '/' '\n'| tail -n 1 )
if [ "$fold" != "exa-bayes" ]; then
    echo "must be executed in exa-bayes folder!"
    exit
fi

if [ $IS_APPLE == 0 ]; then
    readlink=readlink
    ccomp=gcc
    cxxcomp=g++
    system=linux
    compi=mpicc.openmpi
    cxxompi=mpicxx.openmpi
    cmpich=mpicc.mpich2
    cxxmpich=mpicxx.mpich2
else
    export PATH=/opt/local/libexec/gnubin/:/opt/local/bin:/opt/local/sbin:/opt/local/bin:/opt/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin
    readlink=greadlink
    ccomp=clang
    cxxcomp=clang++
    system=apple
    compi=mpicc-openmpi-mp
    cxxompi=mpicxx-openmpi-mp
    cmpich=mpicc-mpich-mp
    cxxmpich=mpicxx-mpich-mp

    export CXXFLAGS="-stdlib=libc++ $CXXFLAGS"
fi




./configure || exit

rm -f  exabayes-*.zip exabayes-*.tar.gz

for vect in $(echo  $ssetypes ) 
do 
    mpis="mpich2 openmpi"
    
    for mpi in $( echo  $mpis )
    do 
	make distclean

	if [ "$mpi" == "mpich2" ]; then
	    mpicxx="$cxxmpich"
	else

	    mpicxx="$cxxompi"
	fi

	arg=""
	if [ "$vect" == "sse" ]; then
	    arg="--disable-avx"
	elif [ "$vect" == "no-sse" ]; then 
	    arg="--disable-sse"
	fi

	mkdir -p  bin && rm -rf bin/*
	mkdir -p  distro-build  && rm -rf distro-build/*

	cd distro-build 

	../configure --enable-mpi  --prefix $($readlink -f ../) CXXFLAGS="-static-libstdc++" CC="ccache $ccomp" CXX="ccache $cxxcomp" MPICXX=$mpicxx $arg ||   exit
	make -j $cores || exit 
	make install  || exit 

	# rm -rf * 

	# ../configure --prefix $($readlink -f ../) LDFLAGS="-static" CC="ccache $ccomp" CXX="ccache $cxxcomp"  $arg ||   exit

	# make -j $cores || exit 
	# make install || exit   

	cd .. 

	# build the distribution 
	./configure 
	make dist 
	make dist-zip 

	# mv packages 
	name=$system-$mpi-$vect
	newname=$(ls exabayes-*.tar.gz | sed "s/\(.*\)\(.tar.gz\)/\1-$name\2/") 
	\mv exabayes*.tar.gz ./packages/$newname
	newname=$(ls exabayes-*.zip | sed "s/\(.*\)\(.zip\)/\1-$name\2/")
	\mv exabayes*.zip ./packages/$newname
    done 
done 

rm bin/*


if [  "$IS_APPLE" == "0" ]; then
    ./configure 
    make dist 
    make dist-zip 
    name=src
    newname=$(ls exabayes-*.tar.gz | sed "s/\(.*\)\(.tar.gz\)/\1-$name\2/")
    \mv exabayes*.tar.gz ./packages/$newname
    newname=$(ls exabayes-*.zip | sed "s/\(.*\)\(.zip\)/\1-$name\2/")
    \mv exabayes*.zip ./packages/$newname
fi
