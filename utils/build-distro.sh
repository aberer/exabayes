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
    compi=mpicc.openmpi
    cxxompi=mpicxx.openmpi
    cmpich=mpicc.mpich2
    cxxmpich=mpicxx.mpich2
else
    compi=mpicc-openmpi-mp
    cxxompi=mpicxx-openmpi-mp
    cmpich=mpicc-mpich-mp
    cxxmpich=mpicxx-mpich-mp
fi

ccomp=gcc
cxxcomp=g++

export OMPI_MPICC=$ccomp
export OMPI_MPIXX=$cxxcomp
export MPICH_CC=$ccomp
export MPICH_CXX=$cxxcomp

./configure 

# make man 
# make distclean

ldflags=" -static-libstdc++"  # -static-libgcc
# ldflags=" -static-libgcc  "

export OMPI_LDFLAGS=$ldflags
export MPICH_LDFLAGS=$ldflags

rm -f  exabayes-*.zip exabayes-*.tar.gz

for vect in $(echo  $ssetypes ) 
do 

    if [ "$vect" == "no-sse" ]
    then 
	mpis=sequential
    else
	mpis="mpich openmpi" 
    fi 

    for mpi in $( echo  $mpis )
    do 
	make distclean

	if [ $mpi == mpich ]; then
	    CC=$cmpich
	    CXX=$cxxmpich
	else
	    CC=$compi
	    CXX=$cxxompi
	fi

	arg=""
	if [ "$vect" == "sse" ]; then
	    arg="--disable-avx"
	elif [ "$vect" == "no-sse" ]; then 
	    arg="--disable-sse"
	fi

	mkdir -p  bin && rm -rf bin/*
	mkdir -p  distro-build  && rm -rf distro-build/*

	if [ "$mpi" != "seqential" ] ; then 
	    # build MPI
	    cd distro-build 
	    # ../configure --enable-mpi  CC="ccache $CC" CXX="ccache $CXX" $arg  --prefix $(readlink -f ../bin) --bindir $(readlink -f  ../bin ) && make -j $cores  && make install
	    cmd="../configure --enable-mpi  CC=\"ccache $CC\" CXX=\"ccache $CXX\" $arg  --prefix $(readlink -f ../bin) --bindir $(readlink -f  ../bin ) && make -j $cores  && make install " 
	    eval $cmd

	    if [ $? != 0 ]; then
		echo -e  "\n\nPROBLEM\n\n"
		echo "$cmd"
		exit
	    fi

	    cd .. 
	fi 

	# build regular 
	rm -rf distro-build/*
	cd distro-build
	../configure  CC="ccache $ccomp" CXX="ccache $cxxcomp"    LDFLAGS="$ldflags" $arg   --prefix $(readlink -f  ../bin) --bindir $(readlink -f ../bin) && make -j$cores   && make install 
	if [ $? != 0  ]; then
	    echo -e "\n\nPROBLEM\n\n"
	    exit
	fi
	cd .. 
	
	# build the distribution 
	./configure 
	make dist 
	make dist-zip 

	# mv packages 
	name=$mpi-$vect
	newname=$(ls exabayes-*.tar.gz | sed "s/\(.*\)\(.tar.gz\)/\1-$name\2/")
	\mv exabayes*.tar.gz ./packages/$newname
	newname=$(ls exabayes-*.zip | sed "s/\(.*\)\(.zip\)/\1-$name\2/")
	\mv exabayes*.zip ./packages/$newname
    done 
done 

rm bin/*

./configure 
# make mydist
make dist 
make dist-zip 
name=src
newname=$(ls exabayes-*.tar.gz | sed "s/\(.*\)\(.tar.gz\)/\1-$name\2/")
\mv exabayes*.tar.gz ./packages/$newname
newname=$(ls exabayes-*.zip | sed "s/\(.*\)\(.zip\)/\1-$name\2/")
\mv exabayes*.zip ./packages/$newname

