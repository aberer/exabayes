#! /bin/sh

rm -rf bin 
mkdir bin 

if [ -f Makefile ]; then
    make distclean
fi

# numCores=`grep "cpu cores" /proc/cpuinfo  | head -n 1 | cut -f 2 -d ':'`

cd bin

# parallel build, if we have additional stuff 

if [  $# -gt 1 ] ; then
    rm -rf build 
    mkdir build
    cd build

    ../../configure --enable-mpi $*

    make

    mv exabayes ..    
    cd .. 

else 
    echo -e "\n\nNo arguments given => will not attempt to build parallel version.\n\n"  
fi 

# cd bin 
rm -rf build 

mkdir build 
cd build 
../../configure 
make

mv yggdrasil postProcParam sdsf credibleSet extractBips consense parser ..

cd .. 
rm -rf build 
