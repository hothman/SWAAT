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


// 		Will  calculate the PSSM for each protein, requires PRODRES
// 		If you use our precalculated PSSM matrices then turn this parameters to False 
// 		In the configuration file (calculate_PSSM = false ) 

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

//		 Alanine scanning  to generate hotspot patch map 
// 		 Will run the alanine scanning prorocol of FoldX over the
//		 PDB structures than will apply an algorithm to locate and 
//		 index the hotspot islands withing the structure 
//		 The PDB file must have the Uniprot ID as basename 
//		 example P04798.pdb
//		 Only proteins defined in the list file 'params.PROTLIST'
//		 Are processed by the workflow


if ( params.calculate_hotspots == true ) {
	// creating output directory 
	// output mapping files to the 'maps' directory
	hotspot_dir  = file("${params.OUTFOLDER}/hotspots")
	hotspot_dir.mkdir() 

	uniprot_list = file("${params.PROTLIST}")
	uniprot_id  = uniprot_list.readLines()
	uniprot_id.remove(0)   // remove the header

	PDBFILES = Channel.fromPath("${params.PDBFILESPATH}/*.pdb") 

	process Hospotislands {
		input:
			val id from uniprot_id
		output: 
			file  '*.csv' 

	    publishDir hotspot_dir , mode:'copy'
	"""
	ln -s ${params.PDBFILESPATH}/${id}.pdb 
	# first repair the structure
	foldx --command=RepairPDB --pdb=${id}.pdb
	# Generate the Ala scan profile
	foldx --command=AlaScan --pdb=${id}_Repair.pdb
	python ${params.SCRIPTHOME}/hotspotPatches.py --pdb ${id}_Repair.pdb --ALAscan *_AS.fxout --suffix ${id}
	"""
	}

}