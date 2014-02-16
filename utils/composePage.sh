#! /bin/sh

if [ $# != 1 ]; then
    echo "content-file"
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
