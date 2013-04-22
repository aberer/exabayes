#! /bin/bash

threshold=90

while true  ; 
do 
    result=$(sudo   sensors | grep Core  | sed 's/[^\+]*+\([0-9]*\)\..*/\1/')

    
    oneTooHot=0
    for elem in $(echo $result ) 
    do	
	if [ "$elem" -gt "$threshold"  ]; then
	    oneTooHot=1
	fi
    done 

    if [ $oneTooHot == 1 ]; then
	killall -s SIGSTOP  exabayes   2> /dev/null
    else
	killall -s SIGCONT  exabayes 2> /dev/null
    fi

    sudo sensors  | grep Core 

    sleep 1     
done 
