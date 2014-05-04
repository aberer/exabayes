#! /bin/bash

topdir=$(dirname  $0 )/../

seed=$RANDOM			# 
# seed=2342
seed=123

# seed=12355

# seed=8106

# src/proposals/
numProc=2
withTree=0

# extraArgs=" -S "
extraArgs=" "  

doParse=0


# early with 150 , VERIFIED 
# seed=31853

# args="--disable-silent-rules" 
# args="--disable-sse"

startFromBest=0
dotests=0

# important: if you do not have google-perftools (and the respective
# *-dev ) package installed, then you should turn this off
useGoogleProfiler=1
useClang=1

cflags=""
cxxflags=""  #  -stdlib=libc++ 

GDB=cgdb


if [ $dotests == 1 ]; then
    args="$args --enable-tests"
fi

# args="$args --disable-silent-rules"
# args="$args" 			#    --disable-sse

runid=testRun

if [ -f /proc/cpuinfo ] ; then 
    numCores=$(cat /proc/cpuinfo  | grep processor  | wc -l) 
else 
    numCores=1

fi 


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
    echo -e  "$0 debug|default|valgrind mpi|thread dataset\n\nwhere the first two arguments are either of the two options, and the third argument is the name of the dataset (e.g., small-dna)"
    exit
fi


mode=$1
codeBase=$2
dataset=$3

# poor...
shift  
shift  
shift  
extra=$*
# echo "extra would be $extra"


pathtodata=$topdir/data/$dataset
if [ ! -d $pathtodata ]; then     
    echo "could not find dataset $dataset"
    exit
fi 


case $mode in
    debug)
	cflags="$cflags -O0 -g"
	cxxflags="$cxxflags -O0 -g"
	gdb="$TERM -e $GDB  -ex run  --args "  #   	
	;;
    default)
	;;
    valgrind)
	cflags="$cflags -O0 -g"
	cxxflags="$cxxflags -O0 -g"
	gdb="$TERM -hold -e valgrind --show-reachable=yes --leak-check=full --tool=memcheck  " #    --leak-check=full --track-origins=yes
	;;
    *)
	echo "mode must be debug, default or valgrind"
	exit
esac


configFile=$pathtodata/config.nex

args="$args --enable-mpi"

if [ "$withTree" = "1" ]; then
    if [ !  -f data/$dataset/tree ]; then
	echo "could not find data/$dataset/tree"
	exit 
    fi

    extraArgs="$extraArgs -t data/$dataset/tree"
    
fi



if [ $doParse -eq 0 ] ; then 
    alnArg=" -f $pathtodata/aln.binary "
else 
    alnArg="-f $pathtodata/aln.phy -q $pathtodata/aln.model"
    if [ ! -f $pathtodata/aln.phy -o ! -f $pathtodata/aln.model  ]; then
	echo "could not find file $pathtodata/aln.{phy,model} "
	exit 
    fi
fi 

if [ "$codeBase" == "mpi" ]; then    
    CC="$ccompiler"  
    CXX="$cxxcompiler"  
    
    baseCall="mpirun -np $numProc  $gdb ./exabayes $alnArg -n $runid -s $seed  $extraArgs -c $configFile $extra"

elif [ "$codeBase" == "thread" ]; then 
    CC="$ccompiler"
    CXX="$cxxcompiler"
    baseCall="$gdb ./yggdrasil -s $seed $alnArg -n $runid $extraArgs -c $configFile $extra "  
else
    echo "second argument must be either 'mpi' or 'thread'"
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
