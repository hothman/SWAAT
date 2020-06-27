#!/usr/bin/python3
__author__ = "Houcemeddine Othman"
__credits__ = "Wits University H3Africa/GSK ADME project"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"

"""
		Scans a directory of PDB files and return which of the chains in the PDBs
		match the input FASTA sequence

			whichPdb.py  --pdbpath path_to_pdbs  --fatsa  sequence.fa --output myoutput.tsv
		
		requires parsePDB.py 

"""

from  parsePDB import ParseFASTA, ParsePDB
from os import path
import glob
import re
import argparse


def process_PDBs(PDBspath, sequence, output="which_pdb.tsv"):
	seq=ParseFASTA(sequence)
	seq.readFASTA()
	seqinFASTA=seq.my_seqs[0]["sequence"]
	if len(seq.my_seqs) >1: 
		raise ValueError("More than one seq in fasta file"  )
	isProt1 = re.match('^[AERTYIPQSDFGHKLMWCVNaertyipqsdfghklmwcvn]+$', seqinFASTA)
	for pdbfile in glob.glob(PDBspath + "*/*.pdb") : 
		mypdb=ParsePDB(pdbfile)
		chains = []
		for chain in mypdb.properties: 
			isProt2=re.match('^[AERTYIPQSDFGHKLMWCVNaertyipqsdfghklmwcvn]+$', chain["seq"] ) 
			if bool(isProt1) and bool(isProt2):
				is_sub_string = seqinFASTA.find(chain["seq"])
				if is_sub_string != -1:
					if chain['chain'] != " ":
						chains.append(chain['chain'])
		#output to file block
		if chains != []:
			with open(output, "a") as outputfile: 
				uniprot_code = seq.my_seqs[0]["header"][0]
				geneID = seq.my_seqs[0]["header"][1]
				basename = path.basename(pdbfile)
				outputfile.writelines("\t".join([geneID, uniprot_code, ",".join(chains), basename+'\n'] )  )


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=" Scans a directory of PDB files and return which of the chains in the PDBs \
		match the input FASTA sequence \n usage: 		whichPdb.py  --pdbpath path_to_pdbs  --fatsa  sequence.fa --output myoutput.tsv")
	parser.add_argument("--pdbpath", help="Path to directory containing PDB files")
	parser.add_argument("--fasta", help="Path sequence file")
	parser.add_argument("--output", help="Output file")
	args = parser.parse_args()

	args.pdbpath != None, 'Wrong usage, try --help'
	args.fasta != None, 'Wrong usage, try --help'

	if args.output == None: 
		process_PDBs(args.pdbpath, args.fasta)
	else: 
		process_PDBs(args.pdbpath, args.fasta, args.output)

