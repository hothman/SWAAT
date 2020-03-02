#!/usr/bin/env nextflow

	/*
	This workflow scans a list of Uniprot IDs in a csv file (params.PROTLIST)
	Download the sequences from ensembl (grch37), the mapping between genomic 
	coordinates and protein coordinates, generate annotation tables, and (optionally)
	calculates PSSMs.

	dependencies :  	fetch_seq.R 		(needs biomaRt and tidyverse)
						prot2genCoor.py 	(needs python3 and transvar installed and configured to use grch37)
						parseProteinAnnotation.py (needs python3)
						You must be connected to www 
	*/


PROTLIST = Channel.fromPath("$params.PROTLIST")
// output sequences to 'sequences' directory 
sequence_dir  = file("${params.OUTFOLDER}/sequences")
sequence_dir.mkdir() 

// fetch sequences 
process GetSequences {
	input: 
		file proteins from PROTLIST
	output:
		file '*.fa' into sequences

	publishDir sequence_dir, mode:'copy'

	"""
	Rscript ${params.SCRIPTHOME}/fetch_seq.R ${params.PROTLIST} ./
	"""
}

/*
sequences.into {seq_for_mapping ; seq_for_annotation }

PROTLIST = Channel.fromPath("$params.PROTLIST")
// output annotation files to the prot_annotation directory
prot_annotation_dir  = file("${params.OUTFOLDER}/prot_annotation")
prot_annotation_dir.mkdir() 

process GetProteinAnnotation {
	input: 
		file(list_of_proteins) from PROTLIST
	output: 
		file '*_annotation.csv' into prot_annotation_file

	publishDir prot_annotation_dir, mode:'copy'

	"""
	tail -n +2  $list_of_proteins >list_of_uniprot_ids
	while read protein_id
	do
		echo \$protein_id
		python ${params.SCRIPTHOME}/parseProteinAnnotation.py --accession \$protein_id

	done < list_of_uniprot_ids
	"""
}

// output mapping files to the 'maps' directory
maps_dir  = file("${params.OUTFOLDER}/maps")
maps_dir.mkdir() 

process GetCoordinates {   

   input:
        file fasta from  seq_for_mapping.flatMap()
   output:
   		file '*.tsv' into data_files	

   	publishDir maps_dir , mode:'copy'
      
   """
   suffix=\$(basename outfolder/${fasta} .fa)
   outputname=\$(echo \$suffix'.tsv') 
   python ${params.SCRIPTHOME}/prot2genCoor.py -i ${fasta} -o \$outputname 
   """
} 

// Will  calculate the PSSM for each protein, requires PRODRES

if ( params.calculate_PSSM == true ) {
    println 'calculate_PSSM'
    // output mapping files to the 'PSSMs' directory
	pssm_dir  = file("${params.OUTFOLDER}/PSSMs")
	pssm_dir.mkdir() 
	// You need psiblast, hmmer and setting PRODRES dependencies 
	process CalculatePSSM {
		input: 
			file fasta from seq_for_annotation.flatMap()
		output:
			file  '*.pssm' into pssm_file

		publishDir pssm_dir , mode:'copy'

		"""
		cat $fasta
		python ${params.PRODRESPATH}/PRODRES.py --input $fasta \
	                  --output ./  \
	                  --pfam-dir ${params.PFAM} \
	                  --pfamscan-script ${params.PRODRESPFAMSCAN} \
	                  --pfamscan_bitscore 2 \
	                  --fallback-db-fasta ${params.UNIREF90} \
	                  --second-search psiblast \
	                  --psiblast_e-val 0.001 \
	                  --psiblast_iter 3 \
	                  --prodres-db ${params.PRODRESDB}

	    rootname=`head -n 1   $fasta |awk -F  "|" '{print \$2}'`
	    outputfilename=\$(echo "\$rootname".pssm )
	    mv \$rootname/outputs/psiPSSM.txt ./\$outputfilename
		"""
	}
}


*/