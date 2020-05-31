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


sequences.into {seq_for_mapping ; seq_for_annotation ; seq_for_chains; seq_for_pdb}

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


/* 		which of the chains in the PDB file correspond to 
		The sequence of the protein.
*/ 
// output mapping files to the 'maps' directory
seq2chain_dir  = file("${params.OUTFOLDER}/Seq2Chain")
seq2chain_dir.mkdir() 

process geneToChainMapping {
	publishDir seq2chain_dir , mode:'copy'
	input:
		file sequence from seq_for_chains.flatMap()
	output: 
		file "${name}_2PDBchain.tsv" into gene2PDBchains
	script: 
		name = sequence.baseName.replaceFirst(".fa","")

	"""
	python  ${params.SCRIPTHOME}/whichPdb.py --pdbpath ${params.PDBFILESPATH}  --fasta ${sequence} --output ${name}_2PDBchain.tsv
	"""
}

/* 		This process generates the mapping between the Uniprot
  		sequence to PDB and calculates sasa, sasa ratio, SS element and H_bonds number
  		per amino acid
*/ 
process uniprot2PDB {
	uniprot2PDB_dir  = file("${params.OUTFOLDER}/uniprot2PDBmap")
	uniprot2PDB_dir.mkdir() 
	publishDir uniprot2PDB_dir , mode:'copy'
	input: 
		file "*" from seq_for_pdb.collect()  //  "*" is used to keep the original names of the files
		file gene2PDBchain from gene2PDBchains
	output: 
		file "*.tsv"

	"""
	gene_name=\$(cut -f 1  $gene2PDBchain)
	pdbfile=\$(cut -f 4  $gene2PDBchain)
	fasta_file=\$(ls \${gene_name}*.fa)
	python ${params.SCRIPTHOME}/parsePDB.py --fasta \$fasta_file  --pdb ${params.PDBFILESPATH}/\$pdbfile
	rm *_2PDBchain.tsv
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

/*		 Alanine scanning  to generate hotspot patch map 
		 Will run the alanine scanning prorocol of FoldX over the
		 PDB structures than will apply an algorithm to locate and 
		 index the hotspot islands withing the structure 
		 The PDB file must have the Uniprot ID as basename 
		 example P04798.pdb
		 Only proteins defined in the list file 'params.PROTLIST'
		 Are processed by the workflow
*/


if ( params.calculate_hotspots == true ) {
	// creating output directory 
	// output mapping files to the 'maps' directory
	hotspot_dir  = file("${params.OUTFOLDER}/hotspots")
	hotspot_dir.mkdir() 

	uniprot_list = file("${params.PROTLIST}")
	uniprot_id  = uniprot_list.readLines()
	uniprot_id.remove(0)   // remove the header


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


/*
		Calculating the covariance matrix and the eigenvectors fors 
		the reference structure using EnCom. 
		The current version of EncoM has an issue with nextflow caused by an exit status
		of 1. To make it work with SWAAT, modify line 370 in src/build_encom.c from "return(1)" with 
		"return(0)" then compile again
*/

PDBLIST = Channel.fromPath("${params.PDBFILESPATH}/*.pdb")

process encomWT {
	publishDir "${params.OUTFOLDER}/ENCoM/", mode:'copy'
	input:
		file pdb from PDBLIST
	output: 
		file "${name}.eige"
		file "${name}.cov"
	script: 
		name = pdb.baseName.replaceFirst(".pdb","")


	"""
	build_encom -i $pdb -cov ${name}.cov -o ${name}.eige 
	"""
}

matrix_dir  = file("${params.OUTFOLDER}/matrices")
matrix_dir.mkdir()

process generate_matrixes {
	publishDir matrix_dir , mode:'copy'

	output:
		file "*.txt"

	"""
	echo "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1" >blosum62.txt

echo "   A R N D C Q E G H I L K M F P S T W Y V
A 0 112 111 126 195 91 107 60 86 94 96 106 84 113 27 99 58 148 112 64
R 112 0 86 96 180 43 54 125 29 97 102 26 91 97 103 110 71 101 77 96
N 111 86 0 23 139 46 42 80 68 149 153 94 142 158 91 46 65 174 143 133
D 126 96 23 0 154 61 45 94 81 168 172 101 160 177 108 65 85 181 160 152
C 195 180 139 154 0 154 170 159 174 198 198 202 196 205 169 112 149 215 194 192
Q 91 43 46 61 154 0 29 87 24 109 113 53 101 116 76 68 42 130 99 96
E 107 54 42 45 170 29 0 98 40 134 138 56 126 140 93 80 65 152 122 121
G 60 125 80 94 159 87 98 0 98 135 138 127 127 153 42 56 59 184 147 109
H 86 29 68 81 174 24 40 98 0 94 99 32 87 100 77 89 47 115 83 84
I 94 97 149 168 198 109 134 135 94 0 5 102 10 21 95 142 89 61 33 29
L 96 102 153 172 198 113 138 138 99 5 0 107 15 22 98 145 92 61 36 32
K 106 26 94 101 202 53 56 127 32 102 107 0 95 102 103 121 78 110 85 97
M 84 91 142 160 196 101 126 127 87 10 15 95 0 28 87 135 81 67 36 21
F 113 97 158 177 205 116 140 153 100 21 22 102 28 0 114 155 103 40 22 50
P 27 103 91 108 169 76 93 42 77 95 98 103 87 114 0 74 38 147 110 68
S 99 110 46 65 112 68 80 56 89 142 145 121 135 155 74 0 58 177 144 124
T 58 71 65 85 149 42 65 59 47 89 92 78 81 103 38 58 0 128 92 69
W 148 101 174 181 215 130 152 184 115 61 61 110 67 40 147 177 128 0 37 88
Y 112 77 143 160 194 99 122 147 83 33 36 85 36 22 110 144 92 37 0 55
V 64 96 133 152 192 96 121 109 84 29 32 97 21 50 68 124 69 88 55 0" >grantham.txt

echo " L I V G A P Q N M T S C E D K R Y F W H
L 0 5 9 24 15 23 22 20 20 23 23 24 30 25 23 33 30 19 30 25
I 5 0 7 25 17 24 24 23 22 21 25 26 31 28 24 34 34 22 34 28
V 9 7 0 19 12 20 25 23 23 17 20 21 31 28 26 36 36 26 37 31
G 24 25 19 0 9 17 32 26 34 20 19 21 37 33 31 43 36 29 39 34
A 15 17 12 9 0 16 26 25 25 20 16 13 34 30 26 37 34 26 36 29
P 23 24 20 17 16 0 33 31 31 25 24 25 43 40 31 43 37 27 37 36
Q 22 24 25 32 26 33 0 10 13 24 21 22 14 22 21 23 29 24 31 27
N 20 23 23 26 25 31 10 0 21 19 15 19 19 14 27 31 28 24 32 24
M 20 22 23 34 25 31 13 21 0 25 22 17 26 31 24 28 32 24 31 30
T 23 21 17 20 20 25 24 19 25 0 12 19 34 29 34 38 32 28 38 34
S 23 25 20 19 16 24 21 15 22 12 0 13 29 25 31 37 29 25 35 28
C 24 26 21 21 13 25 22 19 17 19 13 0 33 28 32 36 34 29 37 31
E 30 31 31 37 34 43 14 19 26 34 29 33 0 7 26 31 34 35 43 27
D 25 28 28 33 30 40 22 14 31 29 25 28 7 0 34 39 34 35 45 35
K 23 24 26 31 26 31 21 27 24 34 31 32 26 34 0 14 34 28 34 27
R 33 34 36 43 37 43 23 31 28 38 37 36 31 39 14 0 36 34 36 31
Y 30 34 36 36 34 37 29 28 32 32 29 34 34 34 34 36 0 13 21 23
F 19 22 26 29 26 27 24 24 24 28 25 29 35 35 28 34 13 0 13 18
W 30 34 37 39 36 37 31 32 31 38 35 37 43 45 34 36 21 13 0 25
H 25 28 31 34 29 36 27 24 30 34 28 31 27 35 27 31 23 18 25 0" >sneath.txt

	"""
}

// stride mutant 