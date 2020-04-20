#!/usr/bin/python3

import re 

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
		source Smitha Dinesh Semwal from www.geeksforgeeks.org"""
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

	def isInProt(self): 
		""" which variants are in the protein coding regions
		The usage of this methods aims to reduce the complexity of 
		the further steps"""
		self.exons()
		print(self.myexons)
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
				 'cDNA_g_codon':cDNA_g_codon ,
				 'AA':AA,
				 'new_codon':new_codon ,
				 'cDNA_new_codon':cDNA_new_codon,
				 'new_AA':new_AA }

	def varToprot(self): 
		self.map.sort(key = lambda x: x[4])

		#print(self.map)
		#sub_li.sort(key = lambda x: x[1])
		for retained_variant in self.retained_vars: 
			#print(self.vars[retained_variant])
			chromosome = self.vars[retained_variant][0]
			position = self.vars[retained_variant][1]
			ref = self.vars[retained_variant][2]
			alt = self.vars[retained_variant][3]
			for var_map in self.map: 
				if  int(var_map[4]) <= position <= int(var_map[5]): 
					position_in_triplet = range( int(var_map[4]) , int(var_map[5])+1 ).index(position)+1
					gDNA_codon = var_map[6] 
					mut_dic = self.mutate(gDNA_codon , alt , position_in_triplet)
					try: 
						var_map[2] == mut_dic['AA']
						print( chromosome, var_map[1], var_map[2] )
						print(mut_dic)
					except:
						pass
					break
			#print(retained_variant)
		#print(self.map)
		#print(self.vars)

		#print(self.getcDNA("GAT"))
		#print(self.mutate('GAT', 'T', 2))



if __name__ == "__main__" : 
	#myvcf = VCF("/home/houcemeddine/BILIM/random_mutation/CYP1A1_dummy.vcf")
	#myvcf.readVcf()

	#mymap = Prot2GenMap("/home/houcemeddine/BILIM/testing_SWAAT/myoutput/maps/CYP1A1.tsv")
	#mymap.readMap()
	myclass = MissenseVars("/home/houcemeddine/BILIM/random_mutation/CYP1A1_dummy.vcf", "/home/houcemeddine/BILIM/testing_SWAAT/myoutput/maps/CYP1A1.tsv"  )
	myclass.isInProt()
	myclass.varToprot()

		