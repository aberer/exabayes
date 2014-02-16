#! /bin/sh

./configure --enable-tests CC="ccache clang" CXX="ccache clang++" && make -j4 && make test 

./configure --enable-tests --enable-mpi CC="ccache mpicc" CXX="ccache clang++"  && make -j4 && make test 


