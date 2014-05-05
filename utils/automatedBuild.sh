#! /bin/bash


apple=administrator@172.16.3.159

./configure 
make distclean

rm packages/* 

# linux build 
rsync --progress -av -C --delete \
    --exclude packages \
    --exclude data \
    --exclude runs \
    --exclude TMP \
    --exclude extra \
     ./ tesla:~/proj/exa-bayes
ssh tesla " cd proj/exa-bayes ; rm packages/*  ;  ./utils/build-distro.sh 0 4 avx sse no-sse ;  "
scp tesla:proj/exa-bayes/packages/* packages/

# apple build 
# rsync --progress -av -C --delete \
#     --exclude packages \
#     --exclude data \
#     --exclude runs \
#     --exclude TMP \
#     --exclude extra \
#     ./ $apple:/Users/administrator/proj/exa-bayes
# ssh $apple " cd proj/exa-bayes ; rm packages/*  ; bash  ./utils/build-distro.sh 1 4 avx sse no-sse ;  "
# scp $apple:/Users/administrator/proj/exa-bayes/packages/* packages/

