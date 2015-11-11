#! /bin/bash


apple=administrator@10.49.49.153

./configure 
make distclean

excludes=" --exclude packages --exclude final-experiments --exclude data --exclude runs --exclude TMP --exclude extra "

# rm packages/* 

# linux build 
# rsync --progress -av -C --delete $excludes ./ tesla:~/proj/exa-bayes
# ssh tesla " cd proj/exa-bayes ; rm packages/* ; ./utils/build-distro.sh 0 4 avx sse no-sse ;  "
# scp tesla:proj/exa-bayes/packages/* packages/

# apple build 

mkdir -p apple-mnt 

rsync --progress -av -C --delete $excludes ./ $apple:/Users/administrator/proj/exa-bayes
ssh $apple " cd proj/exa-bayes ; rm packages/*  ; ./utils/build-distro.sh 1 4 avx sse no-sse ;  "
scp $apple:/Users/administrator/proj/exa-bayes/packages/* packages/
