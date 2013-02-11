#! /bin/bash

if [ "$#" != 3 ]; then
    echo "./plotParamDistribution.sh paramFile pos hist|line" 
    echo "where paramFile is the output of a ExaBayes run and pos is the column in the param file (e.g., 2 for the lnl) "
    exit
fi

theFile=$1
pos=$2
histOpt=$3

if [  "$histOpt"=="hist" -o "$histOpt"=="line" ]; then
    okay=1
else     
    echo "argument 3 must be either hist or line"
fi



tail -n +3 $theFile | cut -f $pos   > tmp 

# echo my dir is  $(dirname $0)
$(dirname $0)/distributionPlotter.R tmp $1.col$2.plot.pdf $histOpt

rm tmp 
