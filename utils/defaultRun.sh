#! /bin/bash


model=GAMMA
seed=4

numCores=$(cat /proc/cpuinfo  | grep processor  | wc -l) 

useClang=1

if [ "$useClang" -ne "0" -a "$(which clang)" != "" ]; then
    ccompiler="clang -Qunused-arguments"
    cxxcompiler="clang++ -Qunused-arguments"
else 
    ccompiler="gcc"
    cxxcompiler="g++"
fi

if [ "$(which ccache)" != "" ]  ; then 
    export CC="ccache $ccompiler"
    export CXX="ccache $cxxcompiler"
else 
    export CC="$ccompiler"
    export CXX="$cxxcompiler"
fi 




if [ "$#" != 3 ]; then
    echo -e  "./defaultRun.sh debug|default pll|examl dataset\n\nwhere the first two arguments are either of the two options, and the third argument is the name of the dataset (e.g., small-dna)"
    exit
fi


args=""

dataset=$3

default=$1
if [ "$default" == "default" ]; then
    args=$args
    gdb=""
elif [ "$default" == "debug" ]; then 
    args="$args --enable-mydebug"
    gdb="$TERM -e gdb -ex run --args "
else  
    echo "first argument must be either 'debug' or 'default'"
    exit 
fi


codeBase=$2
if [ "$codeBase" == "examl" ]; then    
    args="$args --disable-pll"
    baseCall="mpirun -np 2 $gdb ./exabayes -f data/$dataset/aln.examl.binary -n testRun -s $seed -c examples/test.nex" #  -t data/$dataset/tree_2 
elif [ "$codeBase" == "pll" ]; then 
    baseCall="$gdb ./exabayes -s $seed -f data/$dataset/aln.pll.binary -n testRun -c examples/test.nex " #   -t data/$dataset/tree_2 
else
    echo "second argument must be either 'pll' or 'examl'"
    exit
fi

if [ ! -d data/$dataset ]; then
    echo "could not find dataset data/$dataset"
    exit
fi

status="$(./config.status --config | tr -d "'" )"


rm exabayes
./configure -C  $args  $cargs
make -j $numCores

if [ -f ./exabayes ]; then
    echo "calling exabayes as   $baseCall"
    wait 
    $baseCall    
fi


