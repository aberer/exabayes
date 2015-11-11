#! /bin/sh

if [ $# != 1  ]; then
    echo   "USAGE: $0 seq|para\n\nPlease indicate whether you want to find the parallel or sequential version."
    exit
fi


mydir=`dirname $0`

arg=$1
binary=""

case $arg  in 
    para )  binary=exabayes ;; 
    seq ) binary=yggdrasil ;; 
    * ) echo "accepting only 'para' or 'seq' " >&2 ; exit  1 
esac 


result=""
if [ -d  $mydir/../bin -a -f $mydir/../bin/$binary ]; then
    result=$mydir/../bin/$binary
elif [  -d  $mydir/../bin/bin -a -f $mydir/../bin/bin/$binary ] ; then 
    result=$mydir/../bin/bin/$binary
else 
    echo "Could not find exabayes or yggdrasil. Please execute the build.sh script or download a version of this package that contains the executable binaries."  >&2 
    exit 1 
fi

echo $result


