#! /bin/bash 

# this script will compile 

red='\e[0;31m'
NC='\e[0m' # No Color
green='\e[0;32m'

# configure this 
# see  play.h-its.org:/home/andre/compilers for the setup 

GCC_VERSIONS=(4.8 4.9 5 6)
CLANG_VERSIONS=(3.4 3.5 3.6 3.7 3.8 3.9 4.0)

numcores=4

# check gcc versions
for elem in ${GCC_VERSIONS[*]}
do
    WARNLOG=warnlog-gcc$elem-without-mpi.log

    ./configure CC="ccache gcc-$elem" CXX="ccache g++-$elem" > /dev/null
    make clean > /dev/null
    make -j$numcores  > /dev/null 2> $WARNLOG

    if [ $? != 0 ]; then
        echo -e  "[ ${red} failure ${NC} ] with GCC $elem (no MPI)"
    else
        echo -e  "[ ${green} OK ${NC} ] with GCC $elem (no MPI)"
    fi
done

# check gcc versions with MPI
for elem in ${GCC_VERSIONS[*]}
do
    WARNLOG=warnlog-gcc$elem-mpi.log

    ./configure --enable-mpi MPICXX="mpic++.openmpi" OMPI_CXX="ccache g++-$elem" CC="ccache gcc-$elem" CXX="ccache g++-$elem" > /dev/null
    make clean > /dev/null
    make -j$numcores  > /dev/null 2> $WARNLOG

    if [ $? != 0 ]; then
        echo -e  "[ ${red} failure ${NC} ] with GCC $elem (MPI)"
    else
        echo -e  "[ ${green} OK ${NC} ] with GCC $elem (MPI)"
    fi
done

# check clang versions
for elem in ${CLANG_VERSIONS[*]}
do
    WARNLOG=warnlog-gcc$elem-without-mpi.log

    ./configure CC="ccache clang-$elem" \
                CXXFLAGS="-stdlib=libc++" \
                CXX="ccache clang++-$elem" > /dev/null

    make clean > /dev/null
    make -j$numcores > /dev/null 2> $WARNLOG

    if [ $? != 0 ]; then
        echo -e  "[ ${red} failure ${NC} ] with clang-$elem (no MPI)"
    else
        echo -e  "[ ${green} OK ${NC} ] with clang-$elem (no MPI)"
    fi
done

# check clang versions with MPI
for elem in ${CLANG_VERSIONS[*]}
do
    WARNLOG=warnlog-gcc$elem-mpi.log

    ./configure MPICXX="mpic++.openmpi" OMPI_CXX="ccache clang++-$elem" \
                CC="ccache clang-$elem" \
                CXXFLAGS="-stdlib=libc++" \
                CXX="ccache clang++-$elem"  > /dev/null

    make clean > /dev/null
    make -j$numcores > /dev/null 2> $WARNLOG

    if [ $? != 0 ]; then
        echo -e  "[ ${red} failure ${NC} ] with $elem (MPI)"
    else
        echo -e  "[ ${green} OK ${NC} ] with $elem (MPI)"
    fi
done
