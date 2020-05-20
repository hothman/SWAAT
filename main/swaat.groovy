#!/usr/bin/env nextflow
/*
INPUTS: OUTFOLDER Path to output folder
        listmut Path to variants table (could contain multiple lines per gene)  

		Dependencies: Python 3 (scipy, Biopython)
					  FoldX 
*/ 

// home for vcf files 
params.VCFHOME="/home/houcemeddine/BILIM/SWAAT/exampleInputs/vcf/" 
// HOME for maps tsv files
params.MAPHOME="/home/houcemeddine/BILIM/SWAAT/exampleInputs/maps/" 
// Home of the inhouse scripts
params.SCRIPTS="/home/houcemeddine/BILIM/SWAAT/scripts"
// Directory for all outputs 
params.OUTFOLDER = "/home/houcemeddine/BILIM/testing_SWAAT/swaat_output"

// Path to database HOME
params.DATABASE="/home/houcemeddine/BILIM/testing_SWAAT/myoutput"

// list of gene name to process (one gene name per line)
params.GENELIST="../exampleInputs/gene_list.txt"

// home to script file 
params.SCRIPTHOME = "/home/houcemeddine/BILIM/SWAAT/scripts/"

// home to PDBs 
params.PDBFILESPATH = "/home/houcemeddine/BILIM/testing_SWAAT/PDBs" 

// home to PDBs 
params.PDBFILEFIXED = "/home/houcemeddine/BILIM/testing_SWAAT/PDBs" 

// generate a channel from the gene list and replace "\n"  
genes = Channel.fromPath("$params.GENELIST").splitText()  { it.replaceAll("\n", "") }
		.ifEmpty { error "Cannot find genes in gene list" }


//genes.subscribe { println "File: ${it} => ${it}" }
// generate the list of missens variants from VCF and Map files
process generate_swaat_input {
	publishDir "${params.OUTFOLDER}/$gene", mode:'copy'
	input:
		val(gene) from genes.flatMap()
	output: 
		file "${gene}_var2prot.csv" into var2aa_report
		file "${gene}.swaat" into swaat_input
	script: 
		gene_output_dir = file("${params.OUTFOLDER}/$gene")   // create subdirectories (gene names) in the output directory
		gene_output_dir.mkdir()

	"""
	vcf4gene=\$(ls ${params.VCFHOME}/*$gene*.vcf)
	map4gene=\$(ls ${params.MAPHOME}/*$gene*.tsv)
	# gnerate missense variants report and swaat input 
	python ${params.SCRIPTS}/ParseVCF.py --vcf \$vcf4gene \
										 --map \$map4gene  \
										 --output ${gene}_var2prot.csv \
										 --swaat ${gene}.swaat
	"""
}


process filterVars {
	echo true
	/* execution error is ignored if no variant is in the range of the 
	 exons or, it is not covered by the protein structure. or 
	 no *not_processed.tsv is generated 
	*/
	
	input: 
		file(variants) from swaat_input
	output: 
		file("${name}_retained.tsv") optional true into retained_vars
		file("${name}_not_processed.tsv") optional true into gaps_and_non_processed
		//file("${name}_pointer.tsv") optional true into pointerfile
  		val name into bigreceiver

	publishDir "${params.OUTFOLDER}/$name", mode:'copy'

	script:
		name = variants.baseName.replaceFirst(".swaat","")
	"""
	#!/usr/bin/env python
	import glob
	print("$name")
	with open("$variants", "r") as input: 
		lines=input.readlines() 
	for datafile in glob.glob("${params.DATABASE}/uniprot2PDBmap" + "/*.tsv"):
		# supposes that the variant file contains only one gene 
		gene_name_in_vars = lines[0].split()[0]
		with open(datafile, "r") as mydatafile: 
			datalines=mydatafile.readlines() 
			gene_name_in_map=datalines[0].split()[1]
			if  gene_name_in_map == gene_name_in_vars:
				covered_residues=[line.split()[1] for line in datalines[2:]]
				try: 
					info_line = datalines[0].split()
					fasta = info_line[1]+".fasta"
					pdb = info_line[4]
					with open("${name}_pointer.txt", "w") as pointerfile: 
						pointerfile.writelines(fasta+"\t"+pdb)				
				except: 
					pass
				for var in lines: 
					if (var.split()[2] in covered_residues) and (var.split()[3] != "_") :
						with open("${name}_retained.tsv", "a+") as processed_var: 
							processed_var.writelines(var)

					elif var.split()[3] == "_" : 
						with open("${name}_not_processed.tsv", "a+") as not_processed: 
							not_processed.writelines(var)
					else: 
						with open("${name}_not_processed.tsv", "a+") as not_processed: 
							not_processed.writelines(var)
	"""

}


bigreceiver.into {receiver1 ; receiver2; receiver3 }

/*
This process generates the guiding file (which sequence and which PDB for the mutation)
*/

process generateGuidingFile {
	input: 
		//val variant from vars.flatMap()
		val name from receiver1
	output: 
		file "${name}_whichPDB.tsv" into guidingFile

	"""
	python ${params.SCRIPTHOME}/whichPdb.py --fasta     ${params.DATABASE}/sequences/${name}.fa \
											--pdbpath  ${params.PDBFILESPATH} --output ${name}_whichPDB.tsv
	"""
	}


vars = retained_vars.splitText() { it.replaceAll("\n", "") }
//thepointer = pointerfile.splitText() { it.replaceAll("\n", "") }
// foldx 

vars.into { vars4foldx; vars4encom }

process foldX {
	 input:
	 	val variant from vars4foldx.flatMap()
	 	val name from receiver2
	 output: 
	 	file "${name}*.fxout" into mutation_dif_file 
	 	file "*_1.pdb" into mutant_structure_encom, mutant_structure_freesasa, mutant_structure_automute
	 	file "WT_${name}_repaired.pdb" into wt_PDB_repaired

	"""
	echo $variant >afile.txt
	python  ${params.SCRIPTHOME}/processFoldX.py --var afile.txt --seq2chain ${params.DATABASE}/Seq2Chain --output ${name}
	mutation_suffix=\$(sed 's/;//' individual_list_${name}.txt |sed 's/,//')
	pdbfile=\$(cat ${name}_pdb.txt )
	ln -s ${params.PDBFILEFIXED}/\$pdbfile
	foldx --command=BuildModel --pdb=\$pdbfile --mutant-file=individual_list_${name}.txt
	mv Dif_*.fxout ${name}_\${mutation_suffix}.fxout
	mv  WT_*.pdb  WT_${name}_repaired.pdb
	"""
}

// encom 


process encom {
	input: 
		val(variant) from vars4encom.flatMap()
		file(pdb_mutant) from mutant_structure_encom
		file(mutfile) from mutation_dif_file
	output: 
		file "${var_name}.cov" into covarience
		file "${var_name}.eigen"
		val var_name into var_ID1, var_ID2
	 script: 
	 	var_name = mutfile.baseName.replaceFirst(".fxout","")

	"""
	build_encom -i $pdb_mutant -cov ${var_name}.cov -o ${var_name}.eigen
	"""

}

process freesasa {
	input: 
		file(pdb_mutant) from mutant_structure_freesasa
		val var_name from var_ID1
	"""
	freesasa --format seq -n 200 $pdb_mutant >${var_name}.sasa
	"""

}


process stride {
	input: 
		file(pdb_mutant) from mutant_structure_automute
		file(pdb_wt) from wt_PDB_repaired
		val(var_name) from var_ID2
	output: 
		file "Hbond_${var_name}_mut.dat"
		file "Hbond_${var_name}_wt.dat"
	"""
	stride -f $pdb_mutant  -h >Hbond_${var_name}_mut.dat
	stride -f $pdb_wt  -h >Hbond_${var_name}_wt.dat

	"""


}

// automute 


// stride 

//

