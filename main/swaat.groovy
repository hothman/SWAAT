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
// Home to folder containing Blosum62, Grantham and Sneath matrices
params.MATRICES="/home/houcemeddine/BILIM/testing_SWAAT/myoutput/matrices/"
// Path to the pickle file for Random forest prediction 
params.RFMLM="/home/houcemeddine/BILIM/ADME_PGx/SnpsInPdb/MLmodel/swaat_rf.ML"



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
		file("${name}_retained.tsv") optional true into (retained_vars, retained)
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

retained.collect().println()

bigreceiver.into {receiver1 ; receiver2; receiver3 }


vars = retained_vars.splitText() { it.replaceAll("\n", "") }
//thepointer = pointerfile.splitText() { it.replaceAll("\n", "") }
// foldx 
vars.into { vars4foldx; vars4encom ; vars4suffix ; tata; tata2}



/*
This process generates the guiding file (which sequence and which PDB for the mutation)
*/
process generateGuidingFile {
	echo true
	input: 
		val variant from vars4suffix.flatMap()
	output: 
		file "*.id" into id, idpdb
		file "*_whichPDB.tsv" into guidingFile
	publishDir "/home/houcemeddine/BILIM/SWAAT/main/out" , mode:'copy'
	"""
	id=\$(echo $variant |sed 's/ /-/g')
	echo $variant > \$id.id
	gene=\$(echo $variant| awk {'print \$1'}  )
	filename=\$(echo \$gene.fa )
	echo \$filename >mygene
	python ${params.SCRIPTHOME}/whichPdb.py --fasta     ${params.DATABASE}/sequences/\$filename \
											--pdbpath   ${params.PDBFILESPATH} --output \${id}_whichPDB.tsv

	"""
}
	

tata.flatMap().println()
process foldX {
	 input:
	 	val variant from vars4foldx.flatMap()
	 	file(my_id) from id
	 output: 
	 	file "${the_id}_swaat.csv"  into calculated_parameters

	 script:
		the_id = my_id.baseName.replaceFirst(".id","")

	"""
	echo $variant >seq_coord_${the_id}.txt
	python  ${params.SCRIPTHOME}/processFoldX.py --var seq_coord_${the_id}.txt --seq2chain ${params.DATABASE}/Seq2Chain --output ${the_id}
	mutation_suffix=\$(sed 's/;//' individual_list_${the_id}.txt |sed 's/,//')
	pdbfile=\$(cat ${the_id}_pdb.txt )
	Uniprot=\$(basename \$pdbfile .pdb)
	touch \$Uniprot.pointer 
	ln -s ${params.PDBFILEFIXED}/\$pdbfile
	foldx --command=BuildModel --pdb=\$pdbfile --mutant-file=individual_list_${the_id}.txt
	mv Dif_*.fxout ${the_id}_\${mutation_suffix}_suffixed.fxout
	mv  WT_*.pdb  WT_${the_id}_repaired.pdb
	mv *_1.pdb ${the_id}_mutant.pdb


	build_encom -i ${the_id}_mutant.pdb -cov ${the_id}.cov -o ${the_id}.eigen
	
	freesasa -n 200 ${the_id}_mutant.pdb >${the_id}_mut.sasa
	freesasa -n 200 WT_${the_id}_repaired.pdb >${the_id}_wt.sasa
	freesasa --format seq -n 200 WT_${the_id}_repaired.pdb >${the_id}_perAA.sasa 

	stride -f ${the_id}_mutant.pdb -h >Hbond_${the_id}_mut.dat
	stride -f WT_${the_id}_repaired.pdb -h >Hbond_${the_id}_wt.dat

	gene_names=\$(awk {'print \$1'} $my_id)
	coor_in_ref=\$(awk {'print \$1'} seq_coord_${the_id}.txt) 

	ID=\$(basename *.pointer .pointer )
	eigefile=\$(echo \$ID.eige)


	python  ${params.SCRIPTHOME}/parseOutput.py --diff ${the_id}*.fxout  \
											    --matrix ${params.MATRICES}/blosum62.txt \
											    --stride Hbond_${the_id}_mut.dat \
											    --freesasa ${the_id}_mut.sasa \
												--sneath ${params.MATRICES}/sneath.txt \
												--grantham ${params.MATRICES}/grantham.txt \
											    --freesasaWT ${the_id}_wt.sasa \
											    --aasasa ${the_id}_perAA.sasa \
											    --indiv individual_list_${the_id}.txt \
											    --pdbMut ${the_id}_mutant.pdb \
											    --pdbWT WT_${the_id}_repaired.pdb \
											    --modesMut ${the_id}.eigen \
											    --modesWT ${params.DATABASE}/ENCoM/\$eigefile \
											    --pssm ${params.DATABASE}/PSSMs/\$gene_names.pssm \
											    --genename \$gene_names \
											    --output ${the_id}_swaat.csv

	"""
}

/*
process collect {
	echo true
	input: 
		file(dif_file) from mutation_dif_file_collect
		file(indiv) from indiv_file 
		file(mut_hbond) from hbond_mut
		file(wt_hbond) from hbond_wt
		file(sasa_perAA) from perAA_wt_sasa
		file(mut_sasa) from freesasa_mut
		file(wt_sasa) from freesasa_wt
		file(mod_mut) from eigen_mut
		file(mut_structure) from mutant_structure_collection
		file(wt_structure) from wt_PDB_repaired_collection
		file(guide) from guidingFile
		file(coor_in_seq) from seq_coord 
		file(myWTpdb) from pdbid
		file(pointer) from pointer_ID_uniprot
		file(var_id) from varid
	output: 
		file "${the_id}_swaat.csv"  into calculated_parameters
	

	script:
		the_id = var_id.baseName.replaceFirst(".id","")

	"""
	gene_names=\$(awk {'print \$1'} $var_id)
	coor_in_ref=\$(awk {'print \$1'} $coor_in_seq)
	
	ID=\$(basename $pointer .pointer )
	eigefile=\$(echo \$ID.eige)
	cat  $var_id



	python  ${params.SCRIPTHOME}/parseOutput.py --diff $dif_file  \
											    --matrix ${params.MATRICES}/blosum62.txt \
											    --stride $mut_hbond \
											    --freesasa $mut_sasa \
												--sneath ${params.MATRICES}/sneath.txt \
												--grantham ${params.MATRICES}/grantham.txt \
											    --freesasaWT $wt_sasa \
											    --aasasa $sasa_perAA \
											    --indiv $indiv \
											    --pdbMut $mut_structure \
											    --pdbWT $wt_structure \
											    --modesMut $mod_mut \
											    --modesWT ${params.DATABASE}/ENCoM/\$eigefile \
											    --pssm ${params.DATABASE}/PSSMs/\$gene_names.pssm \
											    --genename \$gene_names \
											    --output ${the_id}_swaat.csv


	"""
}

*/
data_vars = calculated_parameters.collectFile(name: "./swaatall" ,  newLine: false, skip: 1, keepHeader: true)


process predict_var  {
	input: 
		file(data4allAavars) from data_vars

	output: 
		file "predicted_outcomes.csv" into predicted_outcomes

	"""
	python ${params.SCRIPTHOME}/deployML.py --inputdata ${data4allAavars} --modelpickle ${params.RFMLM} --output predicted_outcomes.csv 
	"""


} 


process formatReport {
	input: 
		file vars from var2aa_report.collect()

	"""
	for treeFile in ${vars}
	do
		tail -n +2 \$treeFile >>allVariantsInOneFile.csv
	done

	"""

}

