#! /bin/bash

topdir=$(dirname  $0 )/../

seed=$RANDOM

# seed=11436

# seed=31342  			# problematic on tiny-aa
# seed=4045
# seed=1450
# seed=27159
# seed=28978
# seed=6929

# seed=10115 # aa-test
# seed=19180 #aa-test

# seed=2807 # aa-test-one  <= most research so far 

# seed=13290 # tiny-aa

# seed=24066 # problematic with tiny-aa
# seed=32090 # problematic with  143 (on DNA!)

# src/proposals/
numProc=4

# extraArgs="-R 2 -C 2" 
# extraArgs="-M 3 -S  "
# extraArgs="-M 3"
# extraArgs=" -C 4 "
# extraArgs="-C 2"
# extraArgs="-M 3 "
# extraArgs="-M 3 "
# extraArgs="-R 2 " 
# extraArgs="-C 4"
# extraArgs="-S"
# extraArgs="-S  "
# extraArgs="-m"

# early with 150 , VERIFIED 
# seed=31853



# args="--disable-silent-rules" 
# args="--disable-sse"

startFromBest=0
dotests=0

# important: if you do not have google-perftools (and the respective
# *-dev ) package installed, then you should turn this off
useGoogleProfiler=0
useClang=1

if [ $dotests == 1 ]; then
    args="--enable-tests"
fi


# find additional arguments for the call   
# *=$bak

runid=testRun

if [ -f /proc/cpuinfo ] ; then 
    numCores=$(cat /proc/cpuinfo  | grep processor  | wc -l) 
else 
    numCores=1
fi 


# use cgdb, if available 
GDB=gdb
if [ "$(which cgdb )" != ""   ]; then
    GDB="cgdb"
fi
# GDB=gdb


if [ "$useClang" -ne "0" -a "$(which clang)" != "" ]; then
    ccompiler="clang"
    cxxcompiler="clang++"
    cppflags="-Qunused-arguments "  
else 
    ccompiler="gcc"
    cxxcompiler="g++"
fi

if [ $useGoogleProfiler -eq 1  ]; then
    cppflags="$cppflags -D_USE_GOOGLE_PROFILER "
fi


if [ $useGoogleProfiler -eq 1 ]; then
    libs="-lprofiler "
fi


if [ "$#" -lt 3 ]; then
    echo -e  "$0 debug|default pll|examl dataset\n\nwhere the first two arguments are either of the two options, and the third argument is the name of the dataset (e.g., small-dna)"
    exit
fi


dataset=$3

pathtodata=$topdir/data/$dataset
if [ ! -d $pathtodata ]; then     
    echo "could not find dataset $dataset"
    exit
fi 

cflags=""
cxxflags=""  #  -stdlib=libc++ 

default=$1
if [ "$default" == "debug" ]; then 
    cflags="$cflags -O0 -g"
    cxxflags="$cxxflags -O0 -g"
    gdb="$TERM -e $GDB  -ex run  --args "  #   
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
