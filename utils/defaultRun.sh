#! /bin/bash


topdir=$(dirname  $0 )/../

model=GAMMA
seed=12

numCores=$(cat /proc/cpuinfo  | grep processor  | wc -l) 

# important: if you do not have google-perftools (and the respective
# *-dev ) package installed, then you should turn this off
useGoogleProfiler=0		
useClang=1

if [ "$useClang" -ne "0" -a "$(which clang)" != "" ]; then
    ccompiler="clang"
    cxxcompiler="clang++"
    cppflags="-Qunused-arguments -D__STRICT_ANSI__"
else 
    ccompiler="gcc"
    cxxcompiler="g++"
fi

if [ $useGoogleProfiler -eq 1  ]; then
    cppflags="$cppflags -D_USE_GOOGLE_PROFILER"
fi


if [ $useGoogleProfiler -eq 1 ]; then
    libs="-lprofiler"
fi



if [ "$#" != 3 ]; then
    echo -e  "$0 debug|default pll|examl dataset\n\nwhere the first two arguments are either of the two options, and the third argument is the name of the dataset (e.g., small-dna)"
    exit
fi


# args="--disable-silent-rules" 
args=""
dataset=$3

pathtodata=$topdir/data/$dataset
if [ ! -d $pathtodata ]; then     
    echo "could not find dataset $dataset"
    exit
fi 

cflags="-fno-common"
cxxflags="-fno-common"

default=$1
if [ "$default" == "debug" ]; then 
    cflags="$cflags -O0 -g"
    cxxflags="$cxxflags -O0 -g"
    gdb="$TERM -e gdb -ex run --args "
elif [   "$default" != "debug"   -a   "$default" != "default"   ] ; then 
    echo "first argument must be either 'debug' or 'default'"
    exit 
fi


codeBase=$2
if [ "$codeBase" == "examl" ]; then    
    args="$args --disable-pll"

    CC="mpicc -cc=$ccompiler" 
    CXX="mpicxx -cxx=$cxxcompiler"  
    baseCall="mpirun -np 2  $gdb ./exabayes -f $pathtodata/aln.examl.binary -n testRun -s $seed -c $topdir/examples/test.nex"

    # CC="$ccompiler" 
    # CXX="$cxxcompiler"  
    # baseCall="  $gdb ./exabayes -f $pathtodata/aln.examl.binary -n testRun -s $seed -c $topdir/examples/test.nex"
elif [ "$codeBase" == "pll" ]; then 
    CC="$ccompiler"
    CXX="$cxxcompiler"
    baseCall="$gdb ./exabayes -s $seed -f $pathtodata/aln.pll.binary -n testRun -c $topdir/examples/test.nex " 
else
    echo "second argument must be either 'pll' or 'examl'"
    exit
fi

if [ "$(which ccache)" != "" ]  ; then 
    CC="ccache $CC"
    CXX="ccache $CXX"
fi 

args="$args CC=\""$CC"\" CXX=\""$CXX"\""
if [ "$cflags" != "" ]; then
    args="$args CFLAGS=\""$cflags"\""
fi
if [ "$cxxflags" != "" ]; then
    args="$args CXXFLAGS=\""$cxxflags"\""
fi
if [ "$cppflags" != "" ]; then
    args="$args CPPFLAGS=\""$cppflags"\""
fi
if [ "$libs" != "" ]; then
    args="$args LIBS=\""$libs"\""
fi

rm -f exabayes
if [ -f status ] ; then 
    prevStat=$(cat status)
else    
    prevStat=""
fi 

cmd="$topdir/configure $args"

if [ "$prevStat"  == "$cmd" ]  
then     
    echo "=> no need to reconfigure"
else 
    echo "config before: >"$prevStat"<"
    echo "config    now: >$cmd<"

    rm -f  Makefile     
    eval $cmd

    echo "configuring with $cmd"
    if [ -f Makefile ]; then
	echo "$cmd" > status 	
    fi
fi 

make -j $numCores

if [ -f ./exabayes ]; then
    echo "calling exabayes as   $baseCall"
    wait 
    $baseCall    
fi
