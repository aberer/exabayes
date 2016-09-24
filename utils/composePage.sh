#! /bin/sh

if [ $# != 1 ]; then
    echo "content-file"
    exit
fi

read -p "did you follow the following procedure? C-e C-h l p, C-e C-h h h,  make man? [type yes] " res
if [ "$res" != "yes" ]; then
    echo "then do it!"
    exit
fi

input=$1

start=$(grep -i  -n "<body>" $input  | cut -f 1 -d ':')
end=$(grep -i  -n "</body>" $input  | cut -f 1 -d ':')
len=$((end-start))

echo "" > index.html

cat ex-header.html >> index.html
cat $input | tail -n +$((start+1)) | head -n $((len-1))  >> index.html
cat ex-footer.html >> index.html
