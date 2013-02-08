#! /bin/bash


model=GAMMA
seed=123
numCores=4


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
    args="$args --enable-debug"
    gdb="$TERM -e gdb -ex run --args "
else  
    echo "first argument must be either 'debug' or 'default'"
    exit 
fi


codeBase=$2
if [ "$codeBase" == "examl" ]; then    
    baseCall="mpirun -np 2 $gdb ./exabayes -s data/$dataset/aln.examl.binary -t data/$dataset/tree -n testRun -m $model -p $seed -c examples/test.nex "
elif [ "$codeBase" == "pll" ]; then 
    args="$args --enable-pll"
    baseCall="$gdb ./exabayes_pll -T $numCores -p $seed  -s data/$dataset/aln.pll.binary  -t data/$dataset/tree -n testRun -m $model  -c examples/test.nex "
else
    echo "second argument must be either 'pll' or 'examl'"
    exit
fi


if [ ! -d data/$dataset ]; then
    echo "could not find dataset data/$dataset"
    exit
fi


if [ -f lastConfig -a   "$(cat lastConfig)" == "$args"    ]; then 
    echo "no need to re-configure / re-build"
else 
    echo "calling ./configure $args" 
    ./configure $args  
    make clean     
fi

make -j $numCores
echo $args > lastConfig 

echo "calling exabayes as $baseCall"
wait 
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/ncl $baseCall

