#!/bin/bash -ue
tail -n +2  list.csv >list_of_uniprot_ids
while read protein_id
do
	echo $protein_id
	python /home/houcemeddine/BILIM/SWAAT/scripts//parseProteinAnnotation.py --accession $protein_id

done < list_of_uniprot_ids
