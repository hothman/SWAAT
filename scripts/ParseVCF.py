#!/usr/bin/python3




class VCF():
	def __init__(self, vcf):
		self.vcf = vcf

	def readVcf(self):
		with open(self.vcf, 'r') as vcf_file: 
			lines = vcf_file.readlines()
			for line in lines: 
				if "#" not in line:
					print(line.split())







if __name__ == "__main__" : 

	myvcf = VCF("/home/houcemeddine/BILIM/random_mutation/CYP1A1_dummy.vcf")
	myvcf.readVcf()

		