#! /bin/bash 

# this script will compile 

red='\e[0;31m'
NC='\e[0m' # No Color
green='\e[0;32m'

# configure this 
# see  play.h-its.org:/home/andre/compilers for the setup 

GCC_BIN_DIR=../../gcc/bin/
CLANG_DIR=../../clang/

numcores=48

# check gcc versions 
for elem in $(ls $GCC_BIN_DIR)
do 
    ./configure CC="ccache $GCC_BIN_DIR/$elem/bin/gcc" CXX="ccache $GCC_BIN_DIR/$elem/bin/g++" > /dev/null 
    make clean  > /dev/null
    make -j$numcores  > /dev/null

    if [ $? != 0 ]; then
	echo -e  "[ ${red} failure ${NC} ] with $elem"
    else  
	echo -e  "[ ${green} OK ${NC} ] with $elem"
    fi
done 


# check clang versions  
for elem in $(ls $CLANG_DIR)
do 
    ./configure --enable-mpi \
	CC="ccache $CLANG_DIR/$elem/bin/clang" \
	CXX="ccache $CLANG_DIR/$elem/bin/clang++" \
	CXXFLAGS="-stdlib=libstdc++ -nostdinc++ -Wno-unreachable-code -I$GCC_BIN_DIR/gcc-4.9.0/include/c++/4.9.0 -I$GCC_BIN_DIR/gcc-4.9.0/include/c++/4.9.0/./x86_64-unknown-linux-gnu "  \
	CPPFLAGS="-Qunused-arguments" \
	LDFLAGS="-L$GCC_BIN_DIR/gcc-4.9.0/lib64 "   > /dev/null 
    
    make clean > /dev/null 
    make  -j48  > /dev/null 

    if [ $? != 0 ]; then
	echo -e  "[ ${red} failure ${NC} ] with $elem"
    else  
	echo -e  "[ ${green} OK ${NC} ] with $elem"
    fi
done 



