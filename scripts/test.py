#!/usr/bin/python3

"""
This is the implementation of an algorithm to detect the hot spot patches
based on alanine scanning data. User must provide a PDB file and a Foldx alanine scanning 
output file. 3D Patches are composed of at least 3 residues. 
"""

__author__ = "Houcemeddine Othman"
__credits__ = "Wits University H3Africa/GSK ADME project"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"


import argparse
from Bio.PDB.PDBParser import PDBParser 
from Bio.PDB import NeighborSearch, Selection, NeighborSearch
from Bio.PDB.Selection import unfold_entities


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
			'H2S' : 91.9
		   	}

class alaSCanAnalysis():
	def __init__(self, pdb_file):
		self.pdb_file = pdb_file
		parser = PDBParser( QUIET=True )
		self.structure = parser.get_structure('S',self.pdb_file )
		if len(self.structure) >1: 
			warnings.warn('Structure contains more than a model. Only the first one is used')

	def returnChains(self):
		my_structure = self.structure
		if len(my_structure) >  1 :
			raise ValueError('you have multiple models in the PDB file !')
		for model in self.structure :
			monomers = []
			for chain in model: 
				if chain.id != ' ' :
					monomers.append(chain.id)
		return monomers


mystructure = alaSCanAnalysis("../prepare_data/work/f7/c9ec86de786d7f8eeb518479767a02/P05177.pdb")

print( mystructure.returnChains())