#! /bin/bash

cores=4


if [ $# -lt 1 ]; then
    echo -e  "$0 ssetypes..\n\nwhere ssetypes can be  avx,sse, no-sse"
    exit
fi

ssetypes=$*

fold=$(pwd  | tr '/' '\n'| tail -n 1 )
if [ "$fold" != "exa-bayes" ]; then
    echo "must be executed in exa-bayes folder!"
    exit
fi

ccomp=gcc-4.6
cxxcomp=g++-4.6

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
	    ../configure --enable-mpi  CC="ccache mpicc.$mpi" CXX="ccache mpicxx.$mpi" $arg  --prefix $(realpath ../bin) --bindir $(realpath ../bin)    
	    make -j$cores 
	    make install 
	    cd .. 
	fi 

	# build regular 
	rm -rf distro-build/*
	cd distro-build
	../configure  CC="ccache $ccomp" CXX="ccache $cxxcomp"    LDFLAGS="$ldflags" $arg   --prefix $(realpath ../bin) --bindir $(realpath ../bin)  
	make -j$cores 
	make install 
	cd .. 
	
	# build the distribution 
	./configure 
	make mydist 

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
make mydist
name=src
newname=$(ls exabayes-*.tar.gz | sed "s/\(.*\)\(.tar.gz\)/\1-$name\2/")
\mv exabayes*.tar.gz ./packages/$newname
newname=$(ls exabayes-*.zip | sed "s/\(.*\)\(.zip\)/\1-$name\2/")
\mv exabayes*.zip ./packages/$newname

