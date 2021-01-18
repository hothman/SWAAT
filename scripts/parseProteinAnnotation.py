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
		sequence = []
		transcripts_list = []
		for counter,line  in enumerate(raw_data) : 
			is_feature = bool( re.match(r'^FT', line)  ) 
			is_gene_name = bool( re.match(r'^GN\s+Name=', line)  ) 
			is_uniprot_accession = bool( re.match(r'^AC', line)  ) 
			is_refseq_accession = bool( re.match(r'^DR.+RefSeq;', line)  )
			if is_refseq_accession == True: 
				transcripts_list.append(line)

			is_ensembl_accession = bool( re.match(r'^DR\s+Ensembl;', line)  )
			is_sequence = bool( re.match(r'^\s+', line)  )

			# fill the sequence list (each emtpty line in the uniprot file)
			if is_sequence: 
				sequence.append(re.sub(r"[\s]", '', line))

			# get transcript ID
			try: 
				self.transcript

			except: 	
				if is_refseq_accession:
					for element in line.split(): 
						if "NM_" in element: 
							self.transcript = element.split('.')[0]

			try: self.ensembl_id
			except: 
				if is_ensembl_accession:
					ll = line.split() 
					for element in ll: 
						if 'ENSG' in element: 
							self.ensembl_id=element.replace('.','')
			try: 
				self.uniprot_accession
			except:
				if is_uniprot_accession : 
					self.uniprot_accession = line.split()[1].replace(';', '' ) 
			try:
				self.gene_name
			except:
				if is_gene_name: 
					match = re.search(r'Name=\w+', line) 
					self.gene_name = match.group().replace('Name=', '') 


			if is_feature :
				feature_splitted = re.split( "\s\s+" ,line) 
				for feature in KEYS : 
					feature_in_line = bool( re.match(r'^FT[ \t]+%s[ \t]+' % feature, line)  )   
					#print("sdfdfd", feature_in_line)
					if feature_in_line : 
						if "note="  in raw_data[counter+1]:
							specific_annotation = raw_data[counter+1].split("note=")[-1].replace('"', '')
						else: 
							specific_annotation=""
						anotation_line = re.split( "\s\s+" ,line) 
						#print(anotation_line)
						annotation_type = anotation_line[1]
						res_range = anotation_line[-1].split("..")
						start_residue = int( res_range[0])
						try: 
							end_residue = int( res_range[1])
						except: 
							end_residue = int( res_range[0])
						residue_list_with_common_annotation =   list(range(start_residue,end_residue+1) ) 
						for res in residue_list_with_common_annotation: 
						 	self.lines_to_output.append( [str(res), self.gene_name, self.uniprot_accession , features[feature][1], specific_annotation] )

		# The following block will search for the canonical transcript based on the brackets notation: example [Q16348-1]
		# if not it will keep the first transcript exracted in lines 103-107
		if counter == len(raw_data)-1:     # check if we reached the last line in the file
			if len(transcripts_list ) >0 :  # check if the transcripts_list is not empty
				print(transcripts_list)
				for transcript in transcripts_list:
					if "["+self.uniprot_accession+"-1]" in transcript and "NM_" in transcript and "NM_" in transcript :
						splitted = transcript.split(';') 
						for split in splitted : 
							if "["+self.uniprot_accession+"-1]" in split:
								self.transcript = re.match(r'.+NM_[0123456789]+', split).group().replace(" ","") 
			else: 
				warnings.warn( "The parsed file doesn't include a refseq record" )

		# join the sequence together 
		self.sequence = ''.join(sequence) 
		seq_header = ">{0}|{1}|{2}|{3}\n".format(self.uniprot_accession, self.gene_name, self.ensembl_id, self.transcript  )
		self.sequence = seq_header+self.sequence
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
	fasta_filename = myprot.gene_name+".fa"

	print("print sequence file to {}".format(fasta_filename) )	
	with open( fasta_filename, 'w' ) as fasta_file:
		fasta_file.write(myprot.sequence)
	
