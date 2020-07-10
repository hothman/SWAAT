#!/usr/bin/python3
__author__ = "Houcemeddine Othman"
__credits__ = "Wits University H3Africa/GSK ADME project"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"

"""
This script is intended for internal use of swaat
usage: 
		python processFoldX.py --varfile myvarfile.tsv --2pdbchain gene_2PDBchain.tsv --indivoutput indiv.txt --pdbname pdbfile.txt
"""

import glob
import argparse

def readVar(var_input):
	with open(var_input, 'r') as file: 
		mylines1 = file.readlines()
		if len(mylines1) >1 : 
			raise Exception('Only one variant is requested')
		else: 
			return mylines1[0].split() 

def readMap(map_file):
	"""
	Required to convert from the seq coordinates to Pdb coordinates
	"""
	with open(map_file, "r") as mapfile:
		map_lines = mapfile.readlines()
		gene_from_map_file = map_lines[0].split()[1] 
	ID_in_uniprot = []
	ID_in_PDB = []
	for line in map_lines: 
		ID_in_PDB.append(line.split()[1])
		ID_in_uniprot.append(line.split()[2])
	return gene_from_map_file, ID_in_PDB, ID_in_uniprot


def scanFolder(gene_chain_input, vartable, map_file):
	for data_file in glob.glob(gene_chain_input + "/*2PDBchain.tsv"): 
		with open(data_file, "r") as file : 
			gene = file.readlines()[0].split('\t')
			gene_name_in_map, ids_in_pdb , ids_in_uniprot_seq = readMap(map_file)

			if gene[0] == myvar[0]: 
				indiv_table = []
				for chain in gene[2].split(","): 
					# convert from the seq coordinates to Pdb coordinates of the var
					index_var = ids_in_uniprot_seq.index(myvar[2])
					var_id_in_pdb = ids_in_pdb[index_var]
					indiv_table.append(myvar[1]+chain+var_id_in_pdb+myvar[3])
					print(indiv_table)
				indiv_table[-1] = indiv_table[-1]+";\n"
				print(indiv_table)

				return ','.join(indiv_table), gene[-1].replace("\n",'') 

def outputfile(indiv_expression, outputindiv, outputpdbname):
	with open(outputindiv, "w") as outputfile: 
		outputfile.writelines(indiv_expression[0])

	with open(outputpdbname, "w") as outputfile: 
		outputfile.writelines(indiv_expression[1])


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=" Usage: \
	processFoldX.py --varfile myvarfile.tsv --2pdbchain gene_2PDBchain.tsv --indivoutput indiv.txt --pdbname pdbfile.txt")
	# add long and short argument
	parser.add_argument("--var", help="variant file")
	parser.add_argument("--seq2chain", help="Path to seq2chain")
	parser.add_argument("--map", help="Path to mapping file between PDB and sequence")
	parser.add_argument("--output", help="outputname")
	args = parser.parse_args()

	#myvar = readVar("../main/work/60/d2fe89b553151c0adf3119d0216217/myvar.txt")
	myvar = readVar(args.var)
	# myindiv = scanFolder("/home/houcemeddine/BILIM/testing_SWAAT/myoutput/Seq2Chain/", myvar)
	myindiv = scanFolder(args.seq2chain, myvar, args.map)
	# outputfile(myindiv, "output_4swaat.indiv", "pdbfile.txt")
	outputfile(myindiv, "individual_list_"+args.output+".txt", args.output+"_pdb.txt")