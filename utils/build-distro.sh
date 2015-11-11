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
    ccomp=gcc-4.6
    cxxcomp=g++-4.6
    system=linux
    compi=mpicc.openmpi
    cxxompi=mpicxx.openmpi
    cmpich=mpicc.mpich2
    cxxmpich=mpicxx.mpich2
    # ld="LDFLAGS='-static-libgcc -static-libstdc++\'"
else
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
    case "$vect" in 
	"no-sse" ) mpis=sequential ;; 
	"avx" ) mpis="openmpi mpich2" ;;
	"sse" ) mpis="openmpi mpich2" ;;
	*) echo "error: vect must be no-sse, sse or avx" ; exit 
    esac 


    if [ "$vect" == "no-sse" ]
    then 
	mpis="sequential"
    else
	mpis="openmpi mpich2" 
    fi 

    for mpi in $( echo  $mpis )
    do 
	make distclean

	if [ "$mpi" == "mpich2" ]; then
	    cc="$cmpich"
	    cxx="$cxxmpich"
	else
	    cc="$compi"
	    cxx="$cxxompi"
	fi

	arg=""
	if [ "$vect" == "sse" ]; then
	    arg="--disable-avx"
	elif [ "$vect" == "no-sse" ]; then 
	    arg="--disable-sse"
	fi

	mkdir -p  bin && rm -rf bin/*
	mkdir -p  distro-build  && rm -rf distro-build/*

	if [ "$mpi" != "sequential" ] ; then 
	    # build MPI
	    cd distro-build 
	    ../configure --enable-mpi  --prefix $(readlink -f ../) CC="ccache $cc" CXX="ccache $cxx" $arg   && make -j $cores  && make install || exit   
	    # --bindir $(readlink -f  ../ )
	    # eval $cmd

	    cd .. 
	fi 

	# build regular 
	rm -rf distro-build/*
	cd distro-build
	
	
	if [ "$IS_APPLE" == "0" ] ; then 
	    ../configure --prefix $(readlink -f  ../)  CC="ccache $ccomp" CXX="ccache $cxxcomp" $arg LDFLAGS="-static"    && make -j$cores   && make install || exit  
	    # --bindir $(readlink -f ../bin)
	else 
	    ../configure  CC="ccache $ccomp" CXX="ccache $cxxcomp" $arg --prefix $(readlink -f  ../) && make -j$cores   && make install || exit  
	    # --bindir $(readlink -f ../bin)
	fi 

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
