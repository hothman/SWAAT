#!/usr/bin/python3

"""
Parses a VCF file, generates a file containing all missense 
variants using a map file and a file for SWAAT input
"""

__author__ = "Houcemeddine Othman"
__credits__ = "Wits University H3Africa/GSK ADME project"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"


import re 
import os.path
import argparse

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

class VCF():
	def __init__(self, vcf):
		self.vcf = vcf

	def readVcf(self):
		vars = {}
		var_index= 1
		with open(self.vcf, 'r') as vcf_file: 
			lines = vcf_file.readlines()
			for line in lines: 
				if "#" not in line:
					# chromosome, position, ref allele, alternative allele
					vars[var_index] = [line.split()[0], int(line.split()[1]), line.split()[3], line.split()[4] ]
					var_index += 1
		self.vars = vars

class Prot2GenMap(): 
	def __init__(self, mapfile):
		self.mapfile = mapfile
	
	def readMap(self): 
		with open(self.mapfile, 'r') as myfile: 
			 self.map  = [line.split() for line in myfile.readlines()][1:]
		
	def exons(self): 
		"""
		generates intervals of positions corresponding to exons from
		map file 
		"""
		positions= []
		offset = 10
		for amino_acid in self.map: 
			triplet  = list( range( int(amino_acid[4]), int(amino_acid[5])+1 ) ) 
			positions = positions + triplet

		positions.sort() 
		try:
			positions = positions+  list( range(  int(positions[0]-10), int(positions[0]) )) 
			positions.sort() 
			positions = positions + list( range( int(positions[-1]), int(positions[-1])+10  ) )[1:] 
			positions.sort()
		except: 
			pass	
		self.myexons = list( self._interval_extract(positions) )

	def _interval_extract(self, list): 
		""" converts sorted list of integer into intervals 
		source Smitha Dinesh Semwal from www.geeksforgeeks.org
		for internal call by exons(method)"""
		length = len(list) 
		i = 0
		while (i< length): 
			low = list[i] 
			while i <length-1 and list[i]+1 == list[i + 1]: 
				i += 1
			high = list[i] 
			if (high - low >= 1): 
				yield [low, high] 
			elif (high - low == 1): 
				yield [low, ] 
				yield [high, ] 
			else: 
				yield [low, ] 
			i += 1

class MissenseVars(VCF, Prot2GenMap ):
	def __init__(self, vcf, mapfile):
		VCF.__init__(self, vcf)
		Prot2GenMap.__init__(self, mapfile)
		# generate self.vars and self.map
		self.readVcf()		
		self.readMap()

	def checkCompatibilityVCFMAP(self):
		"""
		Make sure that the chromosomes in the map and vcf files are the same
		"""
		assert self.vars[1][0] == self.map[1][0], "Chromosomes in Map and VCF files are not the same"

	def isInProt(self): 
		""" which variants are in the protein coding regions
		The usage of this methods aims to reduce the complexity of 
		the further steps"""
		self.exons()
		self.retained_vars = []
		for variant in self.vars: 
			var_position = self.vars[variant][1] 
			for interval in self.myexons:			 
				if var_position in  range(interval[0], interval[1]+1 )   : 
					self.retained_vars.append(variant)
					break

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
				 'cDNA_g_codon':cDNA_g_codon,
				 'AA':AA,
				 'new_codon':new_codon ,
				 'cDNA_new_codon':cDNA_new_codon,
				 'new_AA':new_AA }

	def varToprot(self, output_file="output_vars.csv"): 
		"""
		a wrapping method that outputs all the non synonomous variants
		to csv file.
		"""
		self.map.sort(key = lambda x: x[4])
		# self.swaat_vars container of variants for SWAAT input 
		self.swaat_vars = []
		#sub_li.sort(key = lambda x: x[1])
		if self.retained_vars: 
			with open( output_file , 'w') as file:
				file.writelines("gene_name,chromosome,position,ref_allele,alt_allele,ref_AA,AA_position,mutant_AA\n")
			for retained_variant in self.retained_vars: 
				chromosome = self.vars[retained_variant][0]
				position = self.vars[retained_variant][1]
				ref = self.vars[retained_variant][2]
				alt = self.vars[retained_variant][3]
				for var_map in self.map: 
					if  int(var_map[4]) <= position <= int(var_map[5]): 
						position_in_triplet = range( int(var_map[4]) , int(var_map[5])+1 ).index(position)+1
						gDNA_codon = var_map[6]

						mut_dic = self.mutate(gDNA_codon , alt , position_in_triplet)
						if var_map[6] == var_map[7]:  # check if the gDNA and the cDNA have the same sequences, if so, use only the information from the tsv file
							mut_dic["cDNA_new_codon"]=mut_dic["new_codon"]
							mut_dic["cDNA_g_codon"] = var_map[7]
							mut_dic["AA"] = DNA_code[mut_dic["cDNA_g_codon"]]
							mut_dic["new_AA"] = DNA_code[mut_dic["cDNA_new_codon"]]

						if  mut_dic["AA"] != var_map[2] : # quality check of the outputs by matching the translated codon with the amino acid in the map file 
							raise Exception("""The translated codon and the amino acid from the map file do not match \n chromosome: {0} position: {1} """.format(chromosome, position )) 

						if mut_dic['new_AA'] != var_map[2]: 
							# in order : gene name, chromosome, position, reference allele, alternative allele, reference AA, AA position, mutant AA
							# output to a file 
							with open( output_file , 'a') as file:
								file.writelines( ','.join([var_map[1], chromosome, str(position), ref, alt,  var_map[2], str(var_map[3]), mut_dic['new_AA']+"\n"]) )
							# fill the list for SWAAT
							data_for_swaat = (var_map[1],var_map[2], str(var_map[3]), mut_dic['new_AA'])
							if  data_for_swaat not in self.swaat_vars:
								self.swaat_vars.append(data_for_swaat)

	def swaatOutput(self, swaat_input="swaat_input.tsv"): 
		if len(self.swaat_vars ) > 0:
			for variant in self.swaat_vars : 
				with open(swaat_input, "a") as file:
					file.writelines( "\t".join([variant[0], variant[1], variant[2],variant[3]+"\n" ] ) )
		else: 
			print("No variants to report")


if __name__ == "__main__" : 
	parser = argparse.ArgumentParser(description="parsing VCF file and generates SWAAT input")
	# Arguments
	parser.add_argument("--vcf", help="Path to the VCF file")
	parser.add_argument("--map", help="Path to the map file")
	parser.add_argument("--output", help="Output file")
	parser.add_argument("--swaat", help="Output file for SWAAT")

	args = parser.parse_args()
	assert args.vcf != None, 'You need at least a VCF and a map file'
	assert args.map != None, 'You need at least a VCF and a map file'
	# instantize
	myclass = MissenseVars(args.vcf, args.map )
	myclass.checkCompatibilityVCFMAP()
	myclass.isInProt()
	try: 
		myclass.varToprot(output_file=args.output)
		print("		output to file {}".format(args.output) )
	except: 
		myclass.varToprot()
		print("	output to file ./output_vars.csv")
	
	try: 
		myclass.swaatOutput(swaat_input=args.swaat)
		print("	output to file {}".format(args.swaat) )
	except: 
		myclass.swaatOutput()
		print("	output to file ./swaat_input.tsv")
		