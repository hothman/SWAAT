#!/usr/bin/python3
__author__ = "Houcemeddine Othman"
__credits__ = "Wits University H3Africa/GSK ADME project"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"

""" 
Reads Unprot entry text file and convert annotation features to csv format
"""

import argparse
import sys
import re
import warnings
from  urllib import request

# check if you have python 3
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

# based on uniprot documentation 
# https://web.expasy.org/docs/userman.html#FT_line
features = { 'SIGNAL':['signal peptide', '1'], 
			'PROPEP': ['propeptide', '2' ],
			'TRANSIT': ['transit peptide', '3'],
			'CA_BIND': ['rsidue part of a calcium binding region', '4'],
			'ZN_FING': ['rsidue part of a zinc binding regon', '5'],
			'DNA_BIND': ['residue implicated in DNA binding', '6'],
			'NP_BIND':['residue implicated in nucleotide phosphate binding, e.g ATP', '7'],
			'REGION': ['residue associated to functional property', '8'],
			'COILED': ['residue belongs to coiled-coil region', '9'],
			'MOTIF':  ['residue belongs to sequence motif of biological interest', '10'],
			'ACT_SITE': ['residue belongs to the active site', '11'],
			'METAL': ['residue belongs of a metal binding site', '12'],
			'BINDING': ['residue belongs to a binding site', '13'],
			'SITE': ['residue involved in a generic activity', '14'],
			'MOD_RES': ['residue involved in a Posttranslational modification', '15'],
			'LIPID': ['residue involved in covalent binding of a lipid moiety', '16'],
			'CARBOHYD': ['residue involved in the attachment with a glycan group', '17'],
			'DISULFID': ['residue involved in a disulfide bond', '18'],
			'CROSSLNK': ['residue involved with a crosslink bondin gwith another amino acid', '19'],
			
			}

KEYS = list( features.keys() ) 

UNIPROT_ADRESS = 'https://www.uniprot.org/uniprot/'

class ParseUniprotAnnotation(object):
	"""docstring for ParseUniprotAnnotation"""
	def __init__(self, uniprot_accession='', uniprot_file=''):
		self.uniprot_accession = uniprot_accession
		self.uniprot_file = uniprot_file

	def _download_uniprot(self):
		url_adress = UNIPROT_ADRESS+self.uniprot_accession+'.txt'
		print('Downloading entry {} ...'.format(self.uniprot_accession))
		uniprot_file = request.urlopen(url_adress).read().decode("utf-8") 
		uniprot_file = uniprot_file.splitlines()
		return uniprot_file
	
	def _read_uniprot(self): 
		with open(self.uniprot_file) as file:
			uniprot_file = file.readlines()
		return uniprot_file 

	def _extract_features(self): 
		for line in  uniprot_file.splitlines() : 
			print(re.match(r'^FT', line) )

	def parse_uniprot(self):
		if self.uniprot_file != '':
			raw_data = self._read_uniprot()
		elif self.uniprot_accession != '':
			raw_data = self._download_uniprot()
		else: 
			raise ValueError("Provide a Uniprot accession or a path to file.")

		self.lines_to_output = []
		for counter,line  in enumerate(raw_data) : 
			is_feature = bool( re.match(r'^FT', line)  ) 
			is_gene_name = bool( re.match(r'^GN', line)  ) 
			is_uniprot_accession = bool( re.match(r'^AC', line)  ) 
			if is_uniprot_accession : 
				self.uniprot_accession = line.split()[1].replace(';', '' ) 
			if is_gene_name : 
				match = re.search(r'Name=\w+', line) 
				gene_name = match.group().replace('Name=', '') 
			if is_feature :
				feature_splitted = re.split( "\s\s+" ,line) 
				for feature in KEYS : 
					if feature in line: 
						 specific_annotation = raw_data[counter+1].split("note=")[-1].replace('"', '') 
						 anotation_line = re.split( "\s\s+" ,line) 
						 annotation_type = anotation_line[1]
						 res_range = anotation_line[-1].split("..")
						 start_residue = int( res_range[0])
						 try: 
						 	end_residue = int( res_range[1])
						 except: 
						 	end_residue = int( res_range[0])
						 residue_list_with_common_annotation =   list(range(start_residue,end_residue+1) ) 
						 for res in residue_list_with_common_annotation: 
						 	self.lines_to_output.append( [str(res), gene_name, self.uniprot_accession , features[feature][1], specific_annotation] )
		return self.lines_to_output, self.uniprot_accession


parser = argparse.ArgumentParser(description=" A script to pull annotations from uniprot file. Provide either the uniprot accession code a path to uniprot text file.")

# add long and short argument
parser.add_argument("--uniprot_file", help="PDB structure of the reference")
parser.add_argument("--accession", help="Uniprot accession")
parser.add_argument("--output_path", help="Path to folder to write the output file")


args = parser.parse_args()

if __name__ == "__main__":
	if not args.output_path:
		output = './'
	else:
		output = args.output_path
	if args.accession != None: 
		myprot=ParseUniprotAnnotation(uniprot_accession = args.accession)
	elif args.uniprot_file != None: 
		myprot=ParseUniprotAnnotation(uniprot_file = args.uniprot_file)
	output_annotation = myprot.parse_uniprot() 
	
	annotation = output_annotation[0] 
	output_name = output+'/'+output_annotation[1]+'_annotation.csv'
	outputheader = ['residue_id', 'gene_symbol', 'uniprot_accession', 'annotation_tag' , 'note' ]
	print('Output to file  {}'.format(output_name))
	with open( output_name , 'w') as file:
		file.writelines( ','.join( outputheader  )+'\n' )
		for line in annotation: 
			file.writelines( ','.join( line  )+'\n' )
	
