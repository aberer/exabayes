#! /bin/bash

if [ "$#" != "1"  ]; then
    echo "./script exabayes-run-id "
    exit
fi


if [ ! -f lib/raxmlHPC-SSE3 ]; then
    echo "Please create a link as follows: ./lib/raxmlHPC-SSE3"
    exit 
fi

runid=$1
raxml=./lib/raxmlHPC-SSE3

bestLnl=$(for elem in $(ls ExaBayes_parameters.$runid.* )  ; do    cut -f 2 $elem  | tail -n +3 | sort -n -r | head -n 1 ; done  | sort -n -r | head -n 1)

line=$(grep -n -h -- "$bestLnl" ExaBayes_parameters.$runid.* |  cut -f 1 -d ':'  )
run=$(grep -- "$bestLnl" ExaBayes_parameters.$runid.*  | cut -f 1 -d ':' | sed 's/.*\.\([0-9]*\)$/\1/')

lineInTopo=$(($line - 2))
./utils/getTopologies.sh ExaBayes_topologies.$runid.$run | sed -n "${lineInTopo}p" 

