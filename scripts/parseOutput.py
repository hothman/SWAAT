#!/usr/bin/python3
__author__ = "Houcemeddine Othman"
__credits__ = "Wits University H3Africa/GSK ADME project"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"

""" Parse workflow and generate text report 
BLOSUM matrix parsing section is a modified version 
of Jeff Wintersinger's code https://gist.github.com/jwintersinger/1870047 """

import argparse
import sys
import numpy as np 
import re
import warnings
from Bio.PDB.PDBParser import PDBParser 
from Bio.PDB import NeighborSearch, Selection, NeighborSearch
from Bio.PDB.Selection import unfold_entities

# check if you have python 3
if sys.version_info[0] < 3:
    raise Exception("Extracts data from the workflow's output \
    	and generates the report")

class Encom:
	"""
	Collects eigenvectors and eigenvalues from 
	encom output file
	"""
	def __init__(self, eigenfile):
		self.eigenfile = eigenfile
		self.eigenvalues = []
		self.eigenvectors =  []

	def parse_eigen(self):
		with open(self.eigenfile, 'r') as file:
			lines = file.read()
		raw_format = lines.split('\n\n')
		eigenvalues = [] ; eigenvectors = []
		for mode in raw_format:
			if "Eigenvalue:" in mode:
		 		header_and_data = re.compile('.*Eigenvalue:|\n').split(mode )  
		 		clean_list = [elem for elem in header_and_data if elem]
		 		eigenvalues.append( float(clean_list[0]) ) 
		 		current_eigenvector = []
		 		for eigenvector_line in header_and_data[1:]:		
		 			component = eigenvector_line.split('\t')
		 			try : 
		 				current_eigenvector.append( float(component[2]) )
		 			except: 
		 				pass
		 		eigenvectors.append( current_eigenvector )
		self.eigenvalues = eigenvalues
		self.eigenvectors = eigenvectors

class Matrix:
  """
  Given a matrix file, extracts the substitution
  score between two AAs
  """ 
  def __init__(self, matrix_filename):
    self._load_matrix(matrix_filename)

  def _load_matrix(self, matrix_filename):
    with open(matrix_filename) as matrix_file:
      matrix = matrix_file.read()
    lines = matrix.strip().split('\n')

    header = lines.pop(0)
    columns = header.split()
    matrix = {}

    for row in lines:
      entries = row.split()
      row_name = entries.pop(0)
      matrix[row_name] = {}
      if len(entries) != len(columns):
        raise Exception('Improper entry number in row')
      for column_name in columns:
        matrix[row_name][column_name] = entries.pop(0)
    self._matrix = matrix

  def lookup_score(self, a, b):
    a = a.upper()
    b = b.upper()
    if a not in self._matrix or b not in self._matrix[a]:
      raise InvalidPairException('[%s, %s]' % (a, b))
    return self._matrix[a][b]

"""
Reference SASA for extended conformations (ALA-X-ALA) of amino acids, calculated
by freesasa from rsa files. 
"""
SASAref = { 'ALA' : 108.76, 'ASN' : 145.15, 
		   	'ARG' : 237.92, 'ASP' : 142.60, 
		   	'CYS' : 132.21, 'GLN' : 178.79, 
			'GLU' : 174.25, 'GLY' : 81.17,
			'HIS' : 183.12, 'ILE' : 175.82,
			'LEU' : 179.53, 'LYS' : 204.96,
			'MET' : 193.15, 'PHE' : 199.76,
			'PRO' : 137.22, 'SER' : 118.33,
			'THR' : 140.63, 'TRP' : 249.14,
			'TYR' : 214.12, 'VAL' : 152.06,
		   	}

# from descriptor JOND750101  AAindex
Hydrophobicity = { 'ALA' : 0.87, 'ASN' : 0.09, 
		   	'ARG' : 0.85, 'ASP' : 0.66, 
		   	'CYS' : 1.52, 'GLN' : 0.00, 
			'GLU' : 0.67, 'GLY' : 0.10,
			'HIS' : 0.87, 'ILE' : 3.15,
			'LEU' : 2.17, 'LYS' : 1.64,
			'MET' : 1.67, 'PHE' : 2.87,
			'PRO' : 2.77, 'SER' : 0.07,
			'THR' : 0.07, 'TRP' : 3.77,
			'TYR' : 2.67, 'VAL' : 1.87,
		   	}

# from descriptor BIGC670101,  AAindex
volume = { 'ALA' : 52.6, 'ASN' : 75.7, 
		   	'ARG' : 109.1, 'ASP' : 68.4, 
		   	'CYS' : 68.3, 'GLN' : 89.7, 
			'GLU' : 84.7, 'GLY' : 36.3,
			'HIS' : 91.9, 'ILE' : 102.0,
			'LEU' : 102.0, 'LYS' : 105.1,
			'MET' : 97.7, 'PHE' : 113.9,
			'PRO' : 73.6, 'SER' : 54.9,
			'THR' : 71.2, 'TRP' : 135.4,
			'TYR' : 116.2, 'VAL' : 85.1,
		   	}

# three letters to one letter code of amino acids
amino_acids = { 'A':'ALA', 'R':'ARG',
				'N':'ASN', 'D':'ASP',
				'C':'CYS', 'Q':'GLN',
				'E':'GLU', 'G':'GLY',
				'H':'HIS', 
				'L':'LEU', 'K':'LYS', 
				'M':'MET', 'F':'PHE', 
				'P':'PRO', 'S':'SER', 
				'T':'THR', 'W':'TRP', 
				'Y':'TYR', 'V':'VAL',
				'I':'ILE'
			}

class collectSASA:
	"""
	read freesasa output file and collect sasa
	value.
	modes:	'all' 		extract sasa of all the protein
			'per_res'	calculate ratio of exposure for one residue
						expects 'chain' and 'id'
	"""
	def __init__(self, sasa_output):
		self.sasa_output = sasa_output
		self._loadSasaOutput()

	def _loadSasaOutput(self):
		with open(self.sasa_output, 'r') as file:
			lines = file.readlines()
		self.lines = lines

	def get_sasa(self, mode, **kwargs): 
		if mode == 'all':
			for line in self.lines: 
				if "Total   :" in line:
					match = re.findall("\d+\.\d+", line )
					SASA = float( match[0] )
			self.SASA = SASA
			if self.SASA == "":
				warnings.warn("Check the output file of freesasa")
			return self.SASA
		elif mode == 'per_res':
			for line in self.lines: 
				if 'SEQ' in line and( int(line.split()[2])  == kwargs['id'] ) and ( line.split()[1]  == kwargs['chain'] ) : 
					sasa_res = float( line.split()[5] )
					res_type = line.split()[3]
					sasa_ref = SASAref[res_type]
					return round(sasa_res/sasa_ref, 3) 

class CollectAutomute:
	"""
	Parse automute output file and returns
	the tag, the confidence score and SS 
	"""
	def __init__(self, automute_output):
		self.automute_output = automute_output
		self._loadSAutoMuteOutput()
	
	def _loadSAutoMuteOutput(self):
		with open(self.automute_output, 'r') as file:
			lines = file.readlines()
		self.line = lines[1].split()
		if len(self.line) != 12:
			warnings.warn('Check the output file of Automute')
	
	def get_mut_data(self):
		if self.line[3] == "Decreased":
			mut_statue = "D"
		elif self.line[3] == "Increased": 
			mut_statue = "I"
		confidence = float(self.line[4])
		sec_structure = self.line[11]
		if "" in ( mut_statue, confidence, sec_structure ):
			warnings.warn("Check the output file of Automute")
		
		return mut_statue, confidence, sec_structure

class ParsePDB:
	"""handle the PDB files and fill and update the properties"""
	def __init__(self, pdb_file):
		self.pdb_file = pdb_file
		parser = PDBParser( QUIET=True )
		self.structure = parser.get_structure('S',self.pdb_file )
	def pdbAttributes(self, chain='A'):
		for my_chain in  self.structure[0] : 
			if my_chain.id == chain: 
				for index, residue in enumerate(my_chain) : 
					if index == 0: 
						id_start_residue =residue.get_full_id()[3][1]
				id_end_residue =id_start_residue+len(my_chain)-1 
		return id_start_residue, id_end_residue
			


# Default heavy atom names for CHARMM27 force field
donors_and_acceptors = {
	'ARG' : ( ('NE', 'NH1', 'NH2'), () ),
	'ASN' : ( ('ND2'), ('OD1') ),
	'ASP' : ( (), ('OD1', 'OD2') ),
	'CYS' : ( ('SG'), ('OD1') ),
	'CYH' : ( (), ('SG') ),
	'GLN' : ( ('NE2'), ('OE1') ),
	'GLU' : ( (), ('OE1', 'OE2') ),
	'HIS' : ( ('ND1', 'NE2'), ('ND1', 'NE2') ),
	'LYS' : ( ('NZ'), () ),
	'MET' : ( (), ('SD') ),
	'SER' : ( ('OG'), ('OG') ),
	'THR' : ( ('OG1'), ('OG1') ),
	'TRP' : ( ('NE1'), () ),
	'TYR' : ( ('OH'), ('OH') ) }

main_chain_donors = ['N']
main_chain_acceptors = ['O', 'OC1', 'OC2']
donors  = ['NE', 'NH1', 'NH2', 'ND2', 'SG' , 'ND1', 'NE2' , 'NZ' , 'OG' , 'OG1', 'NE1' , 'OH' ]
acceptors = ['OD1', 'OD2', 'SG', 'OE1', 'OE2', 'ND1', 'NE2', 'SD', 'OG', 'OG1', 'OH'] 


class Hbonds(ParsePDB):
	"""docstring for HbondRes"""
	def __init__(self, pdb_file):
		super().__init__(pdb_file)
	
	def collectHbondsNumber(self, foldxIndivFile):
		wt_mut_chain_id = whatMutation(foldxIndivFile)	
		prot_chain = wt_mut_chain_id[3]
		prot_mut_id = int(wt_mut_chain_id[2])
		model = self.structure[0]
		myres = model[prot_chain][prot_mut_id]
		atoms  = Selection.unfold_entities(model, 'A') # A for atoms
		ns = NeighborSearch(atoms)
		resname = myres.get_resname()
		h_bond_as_donor = []
		h_bond_as_acceptor = []
		for atom in myres: 
		    if atom.name in donors: 
		        #atoms  = Selection.unfold_entities(model, 'A')
		        #print(atoms)
		        # cutoff of 3.2 angstroms D-A, strong mostly covalent according to
		        # Jeffrey, George A.; An introduction to hydrogen bonding, Oxford University Press, 1997.
		        close_atoms = ns.search(atom.coord, 3.2)
		        for close_atom  in close_atoms:
		            full_atom_id = close_atom.get_full_id() 
		            if (full_atom_id[3][1] != prot_mut_id) and (full_atom_id[4][0] in acceptors+main_chain_acceptors ):
		                acceptor_atom_id = full_atom_id[3][1]
		                h_bond_as_donor.append(acceptor_atom_id)
		    # do the same thing fro acceptor atoms
		    if atom.name in acceptors:
		        close_atoms = ns.search(atom.coord, 3.2)
		        for close_atom  in close_atoms:
		            full_atom_id = close_atom.get_full_id() 
		            if (full_atom_id[3][1] != prot_mut_id) and (full_atom_id[4][0] in donors+main_chain_donors ):
		                acceptor_atom_id = full_atom_id[3][1]
		                h_bond_as_acceptor.append(acceptor_atom_id)
		return len(h_bond_as_donor)+len(h_bond_as_acceptor)

	def collectSaltBridge(self, position, chain):
		model = self.structure[0]	
		myres = model[chain][position]
		atoms  = Selection.unfold_entities(model, 'A') # A for atoms
		ns = NeighborSearch(atoms)
		if myres.get_resname() not in [ 'ARG', 'LYS', 'ASP', 'GLU' ]: 
			is_sb =  0
		else: 
			for atom in myres: 
				atom_id = atom.get_full_id() 
				if atom_id[4][0] in [ 'NH1', 'NH2','NZ' ]:
					close_atoms = ns.search(atom.coord, 3.2) 
					if 'OE1' or 'OE2' in [atomtype.id for atomtype in close_atoms]: 
						is_sb =  1
						break
					else: 
						is_sb =  0
						break
				elif atom_id[4][0] in [ 'OE1', 'OE2' ]:
					close_atoms = ns.search(atom.coord, 3.2)
					if 'NH1' or 'NH2' or 'NZ' in [atomtype.id for atomtype in close_atoms]: 
						is_sb =  1
						break
					else: 
						is_sb =  0
						break

		return is_sb 

residue_order_in_pssm = { 'A':1, 'R':2, 'N':3, 'D':4, 
						  'C':5, 'Q':6, 'E':7, 'G':8, 'H':9, 'I':10, 'L':11,
						  'K':12, 'M':13, 'F':14, 'P':15, 'S':16, 'T':17, 'W':18, 'Y':19, 'V':20 }

class Pssm:
	""" The classe parses a BLAST matrix """
	def __init__(self, PssmFile):	
		with open(PssmFile, 'r') as file: 
				 self.pssm = [line.split()[:-2] for line in file.readlines() if len(line.split() ) == 44 ]
				 
	def _formatPssm(self): 
		self.log_likelihood  = {}
		self.weighted_observed_percentages = {}
		for line in self.pssm : 
			self.log_likelihood[line[0]] =  line[1:22]
			self.weighted_observed_percentages[line[0]] =  list(line[1])+line[22:]
		return self.log_likelihood, self.weighted_observed_percentages 

	def getscore(self, position, mutation): 
		self._formatPssm()
		residue_to_check = self.log_likelihood[str(position)]
		mutation = mutation.upper()
		mutation_order = residue_order_in_pssm[mutation]
		return  residue_to_check[mutation_order] 


def CollectFoldxEnergy(foldx_dif_file): 
	"""
	Parses foldx diff file and extract dG
	"""
	with open(foldx_dif_file, 'r') as file: 
		lines = file.readlines()
	data = lines[-1].split()
	try: 
		total_energy = float(data[1])
	except: 
		warnings.warn("corrupted or non-existent file foldx diff file")
	return total_energy

def allHbonds(stride_output):
	"""
	Extract H bond numbers from stride output
	"""
	with open(stride_output) as file : 
		for line in file:
			if 'HBT' in line: 
				Hb_number = int( line.split()[1] )
		return Hb_number

def CollectAAclass(amino_acid):
	"""
	return the amino acid class
	according to Albatineh and Razeghifard (2008)
	Clustering Amino Acids Using Maximum Clusters Similarity.
	"""
	aa_classes = {'c1': ["M", "V", "I", "L" ] , 
	'c2' :[ "G", "A", "C", "P", "S", "T" ], 
	'c3' :["F", "W", "Y"], 
	'c4': ["R", "H", "K"], 
	'c5': ["N", "Q", "D", "E"], 
	}
	for aa_class in aa_classes : 
		if amino_acid in aa_classes[aa_class]: 
			return aa_class
			break

def calRMSIP(subspace1, subspace2, n = 50):
	""" Calculate the Root Mean square innerproduct
	as defined by Amedei 1999 for n modes """
	if (isinstance(subspace1, Encom) and isinstance(subspace2, Encom)):
		subspace1.parse_eigen()
		sub1 = subspace1.eigenvectors
		subspace2.parse_eigen()
		sub2 = subspace2.eigenvectors
		sum_outer = 0
		for eigenvector1 in sub1[6:6+n]: 
			for eigenvector2 in sub2[6:6+n]:
				sum_outer = sum_outer + (np.dot(np.array( eigenvector1 ), np.array( eigenvector2)) )**2
		return np.sqrt(np.array(sum_outer) /float(n) ) 
	else: 
		raise TypeError ('Arguments must be instances of class Encom')

def getEntropy(wt_subspace, mut_subspace):
	""" Calculates entropy changing between two 
	mutant and wild type"""
	if (isinstance(wt_subspace, Encom) and isinstance(mut_subspace, Encom)):
		wt_subspace.parse_eigen()
		eigen1 = wt_subspace.eigenvalues
		mut_subspace.parse_eigen()
		eigen2 = mut_subspace.eigenvalues
		fraction_eigenvalues = np.array(eigen2[6:])/np.array(eigen1[6:])     # do this or an overflow error will pop-up
		# deltaS_vib 
		return np.log(np.prod(fraction_eigenvalues)) 

def blosumAaClass(foldxIndivFile, matrixFile): 
	"""
	container for Matrix and CollectAAclass
	"""
	with open(foldxIndivFile, 'r') as file: 
		line =  file.readline()
		wt = line[0]
		mut = line[-3]
		matrix = Matrix(matrixFile)
		classMut = CollectAAclass( mut )
		classWT = CollectAAclass( wt )
		return matrix.lookup_score(wt, mut) , classWT, classMut

def whatMutation( foldxIndivFile):
	with open(foldxIndivFile, 'r') as file:
		line =  list(file.readline() )
	line = list( ''.join(e for e in line if e.isalnum()) ) 
	wt = line[0]
	mut = line[-1]
	chain = line[1]	
	del line[0:2]
	del line[-1]
	id =  ''.join(line) 
	return wt, mut , int(id), chain 

def strideSS(chain, position, stride_output): 
	with open(stride_output) as file: 
		lines = file.readlines()
	for line in lines: 
		if 'ASG' in line:
			splitted = line.split() 
			if chain == splitted[2] and str(position) == splitted[3]: 
				return splitted[5]

def getHydrophobicity(foldxIndivFile):
	""" get amino acid hydrophobicity descriptor """
	indiv = whatMutation( foldxIndivFile) 
	wt = amino_acids[indiv[0]]
	mut = amino_acids[indiv[1]]
	return Hydrophobicity[wt], Hydrophobicity[mut]

def getVolume(foldxIndivFile):
	""" get amino acid volume descriptor """
	indiv = whatMutation( foldxIndivFile) 
	wt = amino_acids[indiv[0]]
	mut = amino_acids[indiv[1]]
	return volume[wt], volume[mut]

parser = argparse.ArgumentParser(description=" A script to extract and calculate data from \
									FoldX, ENcom, Automutate2, stride and freesasa.")

# add long and short argument
parser.add_argument("--modesMut", help="EnCOM eigenvectors file of the mutant structure")
parser.add_argument("--modesWT", help="EnCOM eigenvectors file of the wild type structure")
parser.add_argument("--diff", help="Foldx diff file")
parser.add_argument("--indiv", help="Foldx individual file")
parser.add_argument("--matrix", help="Substitution matrix file")
parser.add_argument("--stride", help="stride output file WT")
parser.add_argument("--freesasa", help="freesasa outputfile")
parser.add_argument("--freesasaWT", help="freesasa outputfile WT")
parser.add_argument("--aasasa", help="per residue sasa file")
parser.add_argument("--pdbMut", help="PDB structure of the mutant")
parser.add_argument("--pdbWT", help="PDB structure of the reference")
parser.add_argument("--pssm", help="Path to BLAST PSSM file")
parser.add_argument("--output", help="outputfile")
parser.add_argument("--grantham", help="Grantham matrix file")


args = parser.parse_args()

if not args.output:  
    output = './'
else:
	output = args.output

if __name__ == "__main__":

	chain = whatMutation(args.indiv)[3]
	position = whatMutation(args.indiv)[2]
	wt_res = position = whatMutation(args.indiv)[0]
	mut_res = position = whatMutation(args.indiv)[1]

	try : 
		energy = str( round(CollectFoldxEnergy(args.diff) ,3) )
	except: 
		energy = ""	


	try:		
		# RMSIP
		mut_encom = Encom( args.modesMut )
		mut_encom.parse_eigen()
		mut_eigenvectors = mut_encom.eigenvectors
		mut_eigenvalues = mut_encom.eigenvalues
		wt_encom = Encom( args.modesWT)
		wt_encom.parse_eigen()
		mut_eigenvectors = wt_encom.eigenvectors
		mut_eigenvalues = wt_encom.eigenvalues
		rmsip = str( round( calRMSIP(mut_encom ,  wt_encom) ,3) )
	except: 
		rmsip = ""

	try:
		# entropy 
		entropy = str( round(getEntropy(mut_encom ,  wt_encom) ,3)) 
	except:
	 	entropy = ""

	try: 
		# Blosum and aa classes
		score_blosum_class = blosumAaClass(args.indiv, args.matrix)
		blusom = str( score_blosum_class[0] )
		classwt = score_blosum_class[1]
		classmut = score_blosum_class[2]
	except: 
		blusom = ""
		classwt = ""
		classmut = ""
	try: 
		# grantham distance uses the same method for blosum matrix
		score_grantham_class = blosumAaClass(args.indiv, args.grantham)
		grantham = str( score_grantham_class[0] )

	except: 
		grantham = ""
	try: 
		# sasa
		my_sasa = collectSASA(args.freesasa)
		sasa = str( my_sasa.get_sasa(mode='all') )
	except: 
		sasa = ""

	try: 
		# sasa
		my_sasa = collectSASA(args.freesasaWT)
		sasaWT = str( my_sasa.get_sasa(mode='all') )
	except: 
		sasaWT = ""

	# hydrogen bond of the variant
	try: 
	 	pdb_mut = Hbonds(args.pdbMut)
	 	hb_mut = str(pdb_mut.collectHbondsNumber(args.indiv) )
	 	
	except: 
	 	hb_mut = ''

	 # hydrogen bond of the wt
	try: 
	 	pdb_wt = Hbonds(args.pdbWT)
	 	hb_wt = str(pdb_wt.collectHbondsNumber(args.indiv) )
	except: 
	 	hb_wt = ''

	try: 
		position = whatMutation(args.indiv)[2]
		per_res_sasa = collectSASA(args.aasasa)
		sasa_ratio = str(per_res_sasa.get_sasa(mode='per_res', chain=chain, id=position) ) 
	except:
		sasa_ratio = '' 

	# stride secondary structure
	try: 
		ss = strideSS(chain , position, args.stride)
	except: 
		ss=''

	#salt bridge mut
	try: 
		is_salt_bridge_mut = str( pdb_mut.collectSaltBridge(position, chain) )
	except: 
		is_salt_bridge_mut = '0'

	#salt bridge wt 
	try: 
		is_salt_bridge_wt = str( pdb_wt.collectSaltBridge(position, chain) )
	except: 
		is_salt_bridge_wt = '0'
	try:
		hydrophobcity_wt_mut = getHydrophobicity(args.indiv)
		hydrophobicityWT = str(hydrophobcity_wt_mut[0])
		hydrophobicityMut = str(hydrophobcity_wt_mut[1])
	except:
		hydrophobicityWT = ''
		hydrophobicityMut = ''
	try:
		volume_wt_mut = getHydrophobicity(args.indiv)
		volumeWT = str(volume_wt_mut[0])
		volumeMut = str(volume_wt_mut[1])
	except:
		volumeWT = ''
		volumeMut = ''
	try : 
		mystructure = ParsePDB(args.pdbWT)
		start_res =  mystructure.pdbAttributes(chain=chain)[0]
		shift = start_res -1 
		my_pssm = Pssm(args.pssm)
		pssm_mut = str(my_pssm.getscore( position-shift, mut_res))
		pssm_wt = str(my_pssm.getscore( position-shift, wt_res) )
		print(my_pssm.getscore( position-shift, mut_res),',',my_pssm.getscore( position-shift, wt_res))
	except: 
		pssm_mut = ''
		pssm_wt = ''

	outputlist = [wt_res, mut_res, position ,chain, 
				 energy, 
				 ss, rmsip, entropy,
				 blusom, classwt, classmut, 
				 sasa, sasaWT, 
				 hb_mut, hb_wt, 
				 is_salt_bridge_mut, is_salt_bridge_wt, 
				 sasa_ratio,
				 hydrophobicityWT, hydrophobicityMut,
				 volumeWT, volumeMut, 
				 pssm_mut, pssm_wt ]

	outputheader = ['wt_res','mut_res', 'position', 'chain',  
					'dG', 
					'SecStruc', 'RMSIP', 'dS', 'subScore',
					 'classWT', 'classMut', 'sasa_mut', 
					 'sasa_wt', 'hb_mut', 'hb_wt', 'sb_mut', 
					 'sb_wt', 'sasa_ratio' , 
					 'hyrophob_WT', 'hyrophob_Mut',
					  'volume_WT', 'volume_Mut', 'pssm_mut', 'pssm_wt']

# create csv file 
	with open( output+'/data.txt' , 'w') as file:
		file.writelines( ','.join( outputheader  )+'\n' )
		file.writelines( ','.join( str(item) for item in outputlist  )+'\n' )
