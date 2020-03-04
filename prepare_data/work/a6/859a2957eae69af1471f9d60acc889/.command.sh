#!/bin/bash -ue
suffix=$(basename outfolder/CYP1A1.fa .fa)
outputname=$(echo $suffix'.tsv') 
python /home/houcemeddine/BILIM/SWAAT/scripts//prot2genCoor.py -i CYP1A1.fa -o $outputname
