#! /bin/bash 

if [ $# != 1 ]; then
    echo "$0 <exp>"
    exit
fi

grep -nH -i "$1"  $(find ./src -name "*.[cht]pp") 
