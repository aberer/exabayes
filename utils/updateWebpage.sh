#! /bin/sh

read -p "did you follow the following procedure? C-e C-h l p, C-e C-h h h,  make man? [type yes] " res 
if [ "$res" != "yes" ]; then
    echo "then do it!"
    exit
fi

path=/home/aberer/proj/exa-bayes
if [ "$(pwd)" != "$path"  ]; then
    echo "can only be executed from $path"
    exit
fi

# update packages, if wanted 
read -p "update package? if yes, specify version: " version 
if [ "$version" != "" ]; then
    rsync --no-p --no-g --chmod=Du=rwx,Dgo=rx --progress  -av packages/* dellsco:/scratch/sco/exelixis-new/material/exabayes/$version
    ssh dellsco chmod o+rx /scratch/sco/exelixis-new/material/exabayes/$version
fi

# compose webpages 
cd webpage 
../utils/composePage.sh content.html
cd ../manual/
../utils/composePage.sh content.html
cd .. 

# check if repository is up to date 
cd /home/aberer/proj/exelixis-web
git pull 
cd $path

rsync -C -L  -av webpage/ /home/aberer/proj/exelixis-web/web/software/exabayes/
