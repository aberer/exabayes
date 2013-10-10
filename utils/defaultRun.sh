#! /bin/bash

topdir=$(dirname  $0 )/../


seed=20051 # parsimony prbolem 
# seed=$RANDOM

# seed=$RANDOM
# seed=123
# examples/

numProc=2
# extraArgs="-Q"
# extraArgs="-M 1"
# extraArgs="-m"

# lakner-27
# seed=28233

# small dna 
seed=5594
# seed=9127

startFromBest=0


# find additional arguments for the call   
# *=$bak

runid=testRun
numCores=$(cat /proc/cpuinfo  | grep processor  | wc -l) 


# use cgdb, if available 
GDB=gdb
if [ "$(which cgdb )" != ""   ]; then
    GDB="gdb"
fi


# important: if you do not have google-perftools (and the respective
# *-dev ) package installed, then you should turn this off
useGoogleProfiler=0
useClang=0

if [ "$useClang" -ne "0" -a "$(which clang)" != "" ]; then
    ccompiler="clang"
    cxxcompiler="clang++"
    cppflags="-Qunused-arguments "
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


if [ "$#" -lt 3 ]; then
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
    gdb="$TERM -e $GDB -ex run --args "
elif [   "$default" != "debug"   -a   "$default" != "default"   ] ; then 
    echo "first argument must be either 'debug' or 'default'"
    exit 
fi

codeBase=$2

# poor...
shift  
shift  
shift  
extra=$*
# echo "extra would be $extra"


configFile=$pathtodata/config.nex

if [ "$codeBase" == "examl" ]; then    
    args="$args --enable-mpi"

    CC="mpicc -cc=$ccompiler" 
    CXX="mpicxx -cxx=$cxxcompiler"  
    baseCall="mpirun -np $numProc  $gdb ./exabayes -f $pathtodata/aln.binary -n $runid -s $seed  $extraArgs -c $configFile $extra"

    # CC="$ccompiler" 
    # CXX="$cxxcompiler"  
    # baseCall="  $gdb ./exabayes -f $pathtodata/aln.examl.binary -n $runid -s $seed -c $topdir/examples/test.nex"
elif [ "$codeBase" == "pll" ]; then 
    CC="$ccompiler"
    CXX="$cxxcompiler"
    baseCall="$gdb ./yggdrasil -s $seed -f $pathtodata/aln.binary -n $runid $extraArgs -c $configFile $extra "  
else
    echo "second argument must be either 'pll' or 'examl'"
    exit
fi


if [ $startFromBest == 1 ]; then 
    if [ ! -f  $pathtodata/best.tre ] ; then 
	echo "tried to start from best tree, but could not find $pathtodata/best.tre"
	echo "if you do not have such a tree, deactivate startFromBest" 
	exit 
    fi 
	
    baseCall="$baseCall -t $pathtodata/best.tre"
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

rm -f exabayes yggdrasil
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

if [ -f ./exabayes -o -f ./yggdrasil ]; then
    echo "calling exabayes as   $baseCall"
    rm -f  ExaBayes_*.${runid}*
    wait 
    $baseCall    
fi
