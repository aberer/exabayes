#! /bin/bash


if [ $# != 1 ]; then
    echo "./getTopologies.sh file"
    exit
fi

file=$1

start=$(grep -n translate   $file | cut -f 1 -d ':')
end=$(grep -n tree $file | cut -f 1 -d ':' | head -n 3 | tail -n 1 )

head -n $(($end-1)) $file | tail -n "+$(($start+1))" | tr -d  ",;"  > tmp 

grep tree $file   | tail -n +3  | sed 's/.*\[&U\] \(.*\)/\1/'  > trees


while read line 
do 
    old=$(echo $line | cut -f 1 -d ' '  )
    new=$(echo $line | cut -f 2 -d ' '  )

    echo "s/\([,(]\)${old}:/\1$new:/"
done < tmp   > sedFile 


for elem in $(cat sedFile)
do 
    cat trees | sed "$elem" > bla 
    mv bla trees 
done 

cat trees 

rm tmp sedFile trees 


 
