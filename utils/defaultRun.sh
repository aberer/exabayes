#! /bin/bash


model=GAMMA
seed=1234

numCores=$(cat /proc/cpuinfo  | grep processor  | wc -l) 

useClang=0

if [ "$useClang" -ne "0" -a "$(which clang)" != "" ]; then
    ccompiler="clang -Qunused-arguments"
    cxxcompiler="clang++ -Qunused-arguments"
else 
    ccompiler="gcc"
    cxxcompiler="g++"
fi


if [ "$#" != 3 ]; then
    echo -e  "./defaultRun.sh debug|default pll|examl dataset\n\nwhere the first two arguments are either of the two options, and the third argument is the name of the dataset (e.g., small-dna)"
    exit
fi


args="--disable-silent-rules " # "" # 

dataset=$3

default=$1
if [ "$default" == "debug" ]; then 
    cflags="CFLAGS=\"-O1 -ggdb\""
    cxxflags="CXXFLAGS=\"-O1 -ggdb\""
    gdb="$TERM -e gdb -ex run --args "
elif [   "$default" != "debug"   -a   "$default" != "default"   ] ; then 
    echo "first argument must be either 'debug' or 'default'"
    exit 
fi

codeBase=$2
if [ "$codeBase" == "examl" ]; then    
    CC="mpicc -cc=$ccompiler" 
    CXX="mpicxx -cxx=$cxxcompiler"  
    args="$args --disable-pll"
    baseCall="mpirun -np 2 $gdb ./exabayes -f data/$dataset/aln.examl.binary -n testRun -s $seed -c examples/test.nex"
elif [ "$codeBase" == "pll" ]; then 
    CC="$ccompiler"
    CXX="$cxxcompiler"
    baseCall="$gdb ./exabayes -s $seed -f data/$dataset/aln.pll.binary -n testRun -c examples/test.nex " 
else
    echo "second argument must be either 'pll' or 'examl'"
    exit
fi

if [ "$(which ccache)" != "" ]  ; then 
    CC="ccache $CC"
    CXX="ccache $CXX"
fi 
args="$args CC=\""$CC"\" CXX=\""$CXX"\" $cflags $cxxflags"

if [ ! -d data/$dataset ]; then
    echo "could not find dataset data/$dataset"
    exit
fi

status="$(./config.status --config | tr -d "'" )"

rm -f exabayes


if [ -f status ] ; then 
    prevStat=$(cat status)
else    
    prevStat=""
fi 

cmd="./configure $args"

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
