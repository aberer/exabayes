#!/bin/bash

if [ $# != 3  ]; then
    echo "$0 <file> <line> <exp>"
    exit 
fi

file=$1 
startLine=$2
name=$3

cols=$(sed -n "${startLine}p"  $file  | tr "\t" "\n"| grep -n -i "$name" | cut -f 1 -d ':' | tr "\n" ",")  

cols=$(echo $cols | rev | cut -c 2- | rev )

cut -f "$cols" $file
 
