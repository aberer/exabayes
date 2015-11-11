#! /bin/bash 

if [ $# -lt 2   ]; then
    echo "$0 relBurnin file[..]"
    exit
fi

relBurnin=$1
shift
treeFiles=$*

headtmp=$(tempfile)


# check tree  num 
refnum=$(grep -c "tree gen" $(echo $treeFiles | tr ' ' '\n' | head -n 1 ))
for file in $treeFiles
do 
    num=$(grep -c "tree gen" $file)
    if [ $refnum != $num ]; then
	echo "danger: different number of trees in files"
	exit
    fi
done  

skiptrees=$(echo "(($refnum *  $relBurnin)  + 0.5)  / 1 " | bc)

# check header 
grep  -h -v "tree gen" $(echo $treeFiles | tr ' ' '\n'  | head -n 1  ) | tail -n +3  > $headtmp
for file in $treeFiles
do 
    tmphere=$(tempfile)
    grep -h -v "tree gen" $file  | tail -n +3 > $tmphere

    out=$(diff -q $headtmp $tmphere)

    diff -q $headtmp $tmphere
    
    if [ "$out" != "" ]; then
	echo "danger! files do not follow same taxon ordering. Aborting."

	cat $headtmp
	echo -e  "\n\nAND\n\n" 
	cat $tmphere
	
	exit
    fi
    
    rm $tmphere
done   

# echo "skipping $skiptrees" 1>&2

grep -h -v "tree gen" $(echo $treeFiles | tr ' ' '\n'  | head -n 1  ) | head -n -1

for file in $treeFiles
do 
    cat $file | grep -h  "tree gen" | tail -n +$(($skiptrees+1))
done  

echo "end;"

rm $headtmp
