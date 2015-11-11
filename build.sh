#! /bin/bash

mkdir -p build/{obj-mpi,obj-seq}

if [ -f Makefile ]; then
    make distclean
fi

numCores=$(grep "cpu cores" /proc/cpuinfo  | head -n 1 | cut -f 2 -d ':')

cd build

if [ $# -gt 1 ]; then
    rm -rf obj-mpi/*
    cd obj-mpi  

    ../../configure --enable-mpi $* 
    make -j $numCores

    mv exabayes ..    

else 
    rm -rf obj-seq/*
    cd obj-seq

    ../../configure
    make -j $numCores 

    mv yggdrasil  postProcParam asdsf credibleSet extractBips consense  parser  .. 

fi
