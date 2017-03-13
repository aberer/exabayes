#! /bin/bash

cd $(dirname $0 )

emacs --batch --kill --eval "(progn \
    (find-file \"manual.org\")\
    (org-html-export-to-html) \
    (org-latex-export-to-latex) \
    (org-latex-export-to-pdf))"

bibtex2html --no-abstract --no-keywords library2.bib

pdflatex manual && bibtex manual && pdflatex manual 

# finally correct the citations in html format 



htmlFile=manual.html

for elem in $( egrep -o  "\\\cite\{[^\}]*\}" $htmlFile ) 
do 
    singleEntries=$(echo "$elem" | sed 's/\\cite{\([^}]*\)}/\1/g'  | tr "," " ")

    replacement="["
    isFirst=0
    for entry in $singleEntries
    do 
    	num=$(grep "name=\"$entry\"" library2.html  | sed "s/.*>\([0-9]\)*<\/a.*/\1/")

	if [ "$isFirst" == 1  ]; then
	    replacement=${replacement}","
	fi
	replacement=$replacement$num
	isFirst=1
    done 

    replacement=${replacement}"]"

    escElem=$(echo $elem | sed 's/\\/\\\\\\/g')

    cmd="cat $htmlFile  | sed \"s/$escElem/<a href='manual.html\\#sec-9'>$replacement<\/a>/g\"  > tmp "
    eval $cmd
    
    \mv tmp $htmlFile
done 
