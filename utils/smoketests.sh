#! /bin/bash

# ascAlnList=(bin 010 really-tiny toy-usefull tiny tiny-parted tiny-aa 140-parted 775-short 041 027 027-parted small-dna small-dna-parted small-dna-superparted 050 064 036 043 029 071 059 052 140 140-randpart 354 p-202 143 044 024)

ascAlnList=(bin)

for dataname in ${ascAlnList[*]}
do
    base=data/$dataname
    alnfile=$base/aln.phy

    if [ -f $base/aln.part ]; then
        partcmd="-q $base/aln.part"
    else
        partcmd="-m $(cat $base/model | awk '{print toupper($0)}')"
    fi

    rand=$RANDOM
    cmdline="-f $alnfile $partcmd -s $rand  -c test-config.nex"

    # testing yggdrasil
    cmd="./yggdrasil $cmdline -T 4  -n yggdrasil.$dataname.$rand"
    eval $cmd

    if [ $? -ne 0 ]; then
        echo "WARNING: problem with $cmd" >> ERRORS
    fi

    # testing exabayes
    cmd="mpirun -np 4 ./exabayes $cmdline -n exabayes.$dataname.$rand"
    eval $cmd

    if [ $? -ne 0 ]; then
        echo "WARNING: problem with $cmd" >> ERRORS
    fi

    # testing exabayes
    cmd="mpirun -np 4 ./exabayes $cmdline -R 2 -C 2 -n exabayes.$dataname-para.$rand"
    eval $cmd

    if [ $? -ne 0 ]; then
        echo "WARNING: problem with $cmd" >> ERRORS
    fi
done

# check restart from checkpoint
