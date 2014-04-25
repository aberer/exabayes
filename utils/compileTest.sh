#! /bin/bash 

red='\e[0;31m'
NC='\e[0m' # No Color
green='\e[0;32m'

for elem in 4.6  4.7 4.8 4.9 
do 
    extra=""

    if [ $elem = 4.6 ]; then
	extra="--disable-avx"
    fi

    ./configure CC="ccache gcc-$elem" CXX="ccache g++-$elem" $extra --enable-mpi  > /dev/null
    make clean  > /dev/null
    make -j3  > /dev/null

    if [ $? != 0 ]; then
	echo -e  "[ ${red} failure ${NC} ] with gcc-$elem"
    else  
	echo -e  "[ ${green} OK ${NC} ] with gcc-$elem"
    fi
done 


