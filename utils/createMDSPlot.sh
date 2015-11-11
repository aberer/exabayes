#! /bin/bash

raxml=raxmlHPC 			# TODO set 

if [ "$#" -lt  "2" ]; then
    echo <<EOF   "$0 id thinning [file ... ] 

* id is just an id for the run 
* thinning is a number: use only every n-th sample for plotting 

* [ file ... ] are a bunch of ExaBayes_topology-files.  Make sure that
  they end properly,.e. if you aborted the run, the file will be
  truncated and you have to remove the last line for having a proper
  file."
EOF
    
    exit
fi

id=$1
shift 

thinning=$1
shift

echo -ne  "" >  numTrees.$id
echo -ne  "" > allTopo.$id

idfile=tmp.$id.idfile
echo -ne  "" > $idfile
for  i in $(seq 1 $#)
do
    file=$1
    shift 

    $(dirname $0)/getTopologies.sh $file |  sed -n "0~${thinning}p" | sed -e '$d'   > topoTMP
    cat topoTMP | wc -l  >> numTrees.$id 
    cat topoTMP >> allTopo.$id
    rm topoTMP
    echo $file | cut -f2,3 -d '.' >> $idfile
done 

rm -f  *.rfDistances
$raxml -f r -n rfDistances  -z allTopo.$id  -m GTRCAT  > /dev/null
cat  RAxML_RF-Distances.rfDistances  | cut -f 1,2,3 -d ' ' | tr -d  ':'  > tmp.$id.rf

$(dirname $0 )/createMDSPlotHelper.R tmp.$id.rf numTrees.$id $id $idfile

rm tmp.$id.rf
rm *.rfDistances
# rm tmp.$id.idfile
