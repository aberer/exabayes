#! /bin/bash

if [ ! -d  lib/phylogenetic-likelihood-library/ ]; then
    echo "please install the PLL into lib (=> need the path lib/phylogenetic-likelihood-library/)"
    exit
fi

if [ ! -f utils/setupPLLLinks.sh ]; then
    echo "please call me from the ExaBayes top-level folder (e.g. /home/aberer/proj/ExaBayes)"
    exit
fi

mkdir src/pll
cd src/pll


for elem in $(ls ../../lib/phylogenetic-likelihood-library/*.c )  
do 
    ln -s $elem $(basename $elem .c)-pll.c
done


for elem in $(ls ../../lib/phylogenetic-likelihood-library/*.h )  
do 
    ln -s $elem 
done


mkdir phylip_parser 

cd phylip_parser


for elem in $(ls ../../../lib/phylogenetic-likelihood-library/phylip_parser/*.h  )  
do
    ln -s  $elem 
done

cd ../../../
