#!/usr/bin/python3
__author__ = "Houcemeddine Othman"
__credits__ = "Wits University H3Africa/GSK ADME project"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"


import sys 
import argparse

def getChain(gene_name, refseq_to_PDB_chain_map):
	with open(refseq_to_PDB_chain_map, "r") as file: 
		data = file.readlines()[0].split()
		gene = data[0]
		chains = data[2].split(",")
		assert gene == gene_name, "Gene nema do not match with the map file"
		print(chains[0])

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=" Get chain identifier")
	parser.add_argument("--gene_name", help="Gene name")
	parser.add_argument("--map", help="Path to the file mapping chain ID to gene ID")
	args = parser.parse_args()
	getChain(args.gene_name, args.map)
