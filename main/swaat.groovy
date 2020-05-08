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
	//errorStrategy 'ignore'
	
	input: 
		file(variants) from swaat_input
	output: 
		file("${name}_retained.tsv") optional true into retained_vars
		file("${name}_not_processed.tsv") optional true into gaps_and_non_processed
		//file("${name}_pointer.tsv") optional true into pointerfile
  		val name into receiver

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


/*
This process generates the guiding file (which sequence and which PDB for the mutation)
*/

process generateGuidingFile {
	input: 
		val variant from vars.flatMap()
		val name from receiver
	output: 
		file "${name}_whichPDB.tsv" into guidingFile

	"""
	python ${params.SCRIPTHOME}/whichPdb.py --fasta     ${params.DATABASE}/sequences/${name}.fa \
											--pdbpath  ${params.PDBFILESPATH} --output   ${name}_whichPDB.tsv
	"""
	}


vars = retained_vars.splitText() { it.replaceAll("\n", "") }
//thepointer = pointerfile.splitText() { it.replaceAll("\n", "") }