#! /bin/bash

if [ "$(which sshpass )" == "" ]; then
    echo "please install sshpass"
    exit
fi

apple=administrator@10.49.49.199
applessh="sshpass -p apple@H1TS ssh $apple"
applescp="sshpass -p apple@H1TS scp "

./configure 
make distclean

excludes=" --exclude packages --exclude final-experiments --exclude data --exclude runs --exclude TMP --exclude extra "

rm packages/* 

# linux build 
rsync --progress -av -C --delete $excludes ./ tesla:~/proj/exa-bayes
ssh tesla " cd proj/exa-bayes ; rm packages/* ; ./utils/build-distro.sh 0 4 avx sse no-sse ;  "
scp tesla:proj/exa-bayes/packages/* packages/

# apple build 
rsync --progress -av -C -e "$applessh" --delete $excludes ./ :/Users/administrator/proj/exa-bayes
$applessh " cd proj/exa-bayes ; rm packages/*  ; ./utils/build-distro.sh 1 4 avx sse no-sse ;  "
$applescp $apple:/Users/administrator/proj/exa-bayes/packages/* packages/
