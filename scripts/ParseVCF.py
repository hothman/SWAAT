#!/usr/bin/python3

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
		""" which variants are in the protein coding regions"""
		self.exons()
		print(self.myexons)
		for variant in self.vars: 
			var_position = self.vars[variant][1] 
			for interval in self.myexons:			 
				if var_position in  range(interval[0], interval[1]+1 )   : 
					print(var_position, range(interval[0], interval[1]+1 ) )
					break
		print(self.myexons)

if __name__ == "__main__" : 
	#myvcf = VCF("/home/houcemeddine/BILIM/random_mutation/CYP1A1_dummy.vcf")
	#myvcf.readVcf()

	#mymap = Prot2GenMap("/home/houcemeddine/BILIM/testing_SWAAT/myoutput/maps/CYP1A1.tsv")
	#mymap.readMap()
	myclass = MissenseVars("/home/houcemeddine/BILIM/random_mutation/CYP1A1_dummy.vcf", "/home/houcemeddine/BILIM/testing_SWAAT/myoutput/maps/CYP1A1.tsv"  )
	myclass.isInProt()

		