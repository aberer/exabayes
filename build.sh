#! /bin/bash

mkdir -p build/{obj-mpi,obj-seq}

if [ -f Makefile ]; then
    make distclean
fi

numCores=$(grep "cpu cores" /proc/cpuinfo  | head -n 1 | cut -f 2 -d ':')

rm -f build/* 
cd build
rm -rf obj-mpi/* obj-seq/* 
cd obj-seq

../../configure
make -j $numCores 

mv yggdrasil  postProcParam asdsf credibleSet extractBips consense  parser  .. 

if [ $# -gt 1 ]; then
    cd ../../
    make distclean
    cd build/obj-mpi  
    ../../configure  --enable-mpi $* 
    make -j $numCores
    mv exabayes ..    
fi
