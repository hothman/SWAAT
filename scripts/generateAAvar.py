#!/usr/bin/python3
#!/usr/bin/python3
__author__ = "Houcemeddine Othman"
__credits__ = "Wits University H3Africa/GSK ADME project"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"

"""
This script scan a vcf file and returns the list of non synonymous
amino acid variants, as input for the workflow. If the variant is in the 
splicing boundaries, it returns a warning and a lsit of these variants. 
Written with python3.5, argument parsing may not work with newer version.
  mode: 
	"report" (default) : only report the non synonymous variants of amino acids ( a tsv file).
	"extended": report a more detailed output (extension *extended.tsv) 
				and the default report

Usage:
generateAAvar.py -i<input-vcf> -d <database-file> [-o <output-file>] [-m <mode>]
"""

import sys
import re
import getopt
import os
import subprocess
import warnings

# check if you have python 3
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

       ###  Check the arguments ###
args, out=getopt.getopt(sys.argv[1:], 'i:o:d:m:' )
argdic={}
for elem in args:
    argdic[elem[0]] = elem[1]

# Check the input vcf
try:
    vcfInput = argdic['-i']
except:
    sys.exit('No VCF file is provided. We stop here.\
    \nUsage:\n\tgenerateAAvar.py -i<input-vcf> -d <database-file> -o <output-file> ')

# Check the database 
try:
    database = argdic['-d']
except:
    sys.exit('Missing database file. We stop here.\
    \nUsage:\n\tgenerateAAvar.py -i<input-vcf> -d <database-file> -o <output-file> ')

# Check the output suffix
try:
    outputsuffix = argdic['-o']
except:
    pass

# Check the mode 
try:
    mode = argdic['-m']
except:
    mode = "report"

print(argdic)
DNA_code = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 

class ParseVCF():
	def __init__(self, vcf_file):
		self.vcf_file = vcf_file

		try:
			myfile = open(self.vcf_file, 'r')
			self.vcf_lines = myfile.readlines()
		finally:
			myfile.close()
		self.variants = []

	def getData(self):
		""" Parse the lines of the VCF file """
		regexp1 = re.compile(r'^#')
		variants = []
		for line in self.vcf_lines[1:]:     # skip the first line
			if regexp1.match(line) is None :
				var_line = line.split('\t')
				chrom = var_line[0]
				POS = var_line[1]
				ID = var_line[2]
				REF = var_line[3]
				ALT = var_line[4]
				variants.append([int(chrom), int(POS), ID, REF, ALT])
		self.variants = variants

	def var2Tsv(self, outputfile):
		"""Convert the VCF vanriant table to TSV"""
		with open( outputfile , 'w') as file:
			file.writelines("chr\tPOS\tID\tREF\tALT\n")
			file.writelines(str(variant[0])+'\t'+ str(variant[1])+'\t' \
			 + variant[2]+'\t'+variant[3]+'\t'+variant[4]+'\n' for variant in self.variants)  

		
class GenerateVarInput():
	"""formats and generates the variants table for the workflow"""
	def __init__(self, vcf_file, database_path):
		self.vcf_file = vcf_file
		self.vcf = ParseVCF(self.vcf_file)

		# read the database file 
		try:
			myfile = open(database_path, 'r')
			self.database = myfile.readlines()
		finally:
			myfile.close()

	def getcDNA(self, DNAseq):
		"""
		returns the cDNA sequence in 5' 3' order 
		"""
		complementary_nuc = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
		isDNA = re.match('^[CGTAcgta]+$', DNAseq)
		compDNA35 = []
		if bool(isDNA )  == True:
			for nucleotide in DNAseq: 
				compDNA35.append( complementary_nuc[nucleotide] )

			compDNA35 =  ''.join(compDNA35)  
			compDNA53 = ''.join(  list(reversed(compDNA35 )  )   ) 
			return  compDNA53 
		else: 
			print('The DNA sequence contains non A, C, T, G characters')


	def mutate(self,g_codon, var, position):
		"""
		returns a dictionary of WT and mutated condons 
		and corresponding amino acids 
			codon: WT codon
			var: mutant nucleotide
			position: 1, 2 or 3
		"""
		global DNA_code
		cDNA_g_codon = self.getcDNA(g_codon)
		AA = DNA_code[cDNA_g_codon] 
		new_codon = list(g_codon)
		new_codon[position-1] = var
		new_codon = ''.join(new_codon)
		cDNA_new_codon = self.getcDNA(new_codon)
		new_AA = DNA_code[cDNA_new_codon] 
		return { 'g_codon':g_codon, 
				 'cDNA_g_codon':cDNA_g_codon ,
				 'AA':AA,
				 'new_codon':new_codon ,
				 'cDNA_new_codon':cDNA_new_codon,
				 'new_AA':new_AA }


	def getVars(self, mode="report", suffix=''):
		global DNA_code
		self.vcf.getData()
		if suffix == '' :
			outputfile = "outputAAvar.tsv"
			fulloutput = "outputAAvar_extended.tsv"
		else: 
			outputfile = suffix+'.tsv'
			fulloutput = suffix+'_extended.tsv'
		with open( outputfile , 'w') as file:
			file.writelines('HGCN'+'\t'+'new_AA'+'\n' )
		if mode =='extended':
			with open( fulloutput , 'w') as file:
				file.writelines('HGCN'+'\t'+
						  		'ref_codon'+'\t'+
						  		'Ref_AA'+'\t'+
						  		'variant_codon'+'\t'+
						  		'variant_AA'+'\t'+
						  		'genome_position'+'\t'+
						  		'\n' )
		for variant in self.vcf.variants:
			alternative_var = variant[4]
			reference_var = variant[3]
			for amino_acid in self.database[1:]:
				my_AA = amino_acid.split('\t')
				chromosome_var = variant[0]
				chromosome_AA = int( my_AA[0] )
				AA_position = int( my_AA[3] )
				hgcn = my_AA[1]
				if chromosome_var  == chromosome_AA  :
					gnom_position = variant[1]
					codon_start = int( my_AA[4] )
					codon_end = int( my_AA[5] )
					# If variant in exon
					if gnom_position in list( range(codon_start, codon_end+1 )  ): 
						triplet = list( range(codon_start, codon_end+1 )  )
						codon_idx={ triplet[0]: 1,
						            triplet[1]: 2, 
						            triplet[2]: 3}
						WT_codon = my_AA[6]   
						var_position_in_codon = codon_idx[gnom_position]
						snv_aa = self.mutate(WT_codon, alternative_var, var_position_in_codon )
						if mode =="report":
							if snv_aa['AA'] != snv_aa['new_AA']:
								with open( outputfile , 'a') as file:
									file.writelines( hgcn+'\t'+snv_aa['new_AA']+str(AA_position)+'\n' )
						if mode =="extended":
							if snv_aa['AA'] != snv_aa['new_AA']:
								with open( outputfile , 'a') as file:
									file.writelines( hgcn+'\t'+snv_aa['new_AA']+str(AA_position)+'\n' )
							with open( fulloutput , 'a') as file:
								#print( snv_aa )
								file.writelines( hgcn+'\t'+
												     snv_aa['g_codon']+'\t'+
													 snv_aa['AA']+'\t'+
													 snv_aa['new_codon']+'\t'+
													 snv_aa['new_AA']+'\t'+
													 str(gnom_position)+'\t'+
													 '\n' )
					elif gnom_position in [codon_start, codon_start+1, codon_end-1, codon_end  ] : 
						with open( fulloutput+".boundaries" , 'a') as file:
							warning_msg = "Nucleotide at position "+str(gnom_position)+" chromosome "+str(chromosome_var )+ \
											"is at a boundary of an exon. \nSee file "+fulloutput+".boundaries"
							file.writelines(  hgcn+'\t'+
											  str(chromosome_var )+'\t'+
											  str(gnom_position)+
											  '\n')


# The workflow 
# Create an instance of GenerateVarInput
processvcf = GenerateVarInput( vcfInput , database ) 
# generate the report 
processvcf.getVars(mode=mode,  suffix=outputsuffix)
