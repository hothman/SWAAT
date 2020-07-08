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
			'HIS' : 91.9, 'ILE' : 102.0, 'H1S' : 91.9,
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

	def _returnChains(self):
		my_structure = self.structure
		if len(my_structure) >  1 :
			raise ValueError('you have multiple models in the PDB file !')
		for model in self.structure :
			monomers = []
			for chain in model: 
				if chain.id != ' ' :
					monomers.append(chain.id)

		return monomers
		
	def readALA(self, ala_scan_file): 
		""" read the alanine scanning file a
		'self.dic' a dictionary where  keys are amino acid positions and
		values are tuples  of amino acid Id and energy """
		self._chainWalk()
		with open(ala_scan_file) as file: 
			lines = file.readlines()
		self.dic = {}
		ala_scan_lines = []
		for idx, line in enumerate(lines): 
			splitted = line.split()
			amino_acid =  splitted[0]
			position = splitted[1] 
			energy = float(splitted[7])
			assert position != self.chainId[idx][0] , 'Residue numbering conflict between PDB (residue {0}) file and alanine scanning file (residue {1})'.format(self.chainId[idx][0],position  )
			self.dic[position] = ( amino_acid , energy, self.chainId[idx][1], self.chainId[idx][0] )
			ala_scan_lines.append( (amino_acid , energy, self.chainId[idx][1], self.chainId[idx][0] ) ) 
		return ala_scan_lines

	def _chainWalk(self): 
		chainId = []
		if len(self.structure) >  1 :
			raise ValueError('you have multiple models in the PDB file!')
		for chain in self.structure[0]: 
			for res in chain: 
				if res.id[0] == ' ':
					chainId.append( (res.id[1], chain.id) )
		self.chainId = chainId

	def allChainsClusters(self, ala_scan_file):
		alascan_table = self.readALA(ala_scan_file)	
		All_chains_dics = {}
		for chain in self._returnChains():
			chain_dic = {}
			for line in alascan_table: 
				if chain == line[2]: 
					position = line[3]
					amino_acid = line[0]
					energy = line[1]
					chain_dic[position] = ( amino_acid , energy, chain, position )
			All_chains_dics[chain] = chain_dic
		return  All_chains_dics

	def DefinePatches(self, chain, dic_chain,energy_cutoff=2.0,  dist_cutoff=6):
		"""
		This method is the implementation of the algorithm
		"""
		self.dic = dic_chain 
		clusters_pre_list = []     # The list containing the clusters
		self.chain = chain
		if 'dic'  in dir(self):
			self.protein_ca_atoms = [atom for atom in self.structure[0].get_atoms() if atom.name=="CA"]
			self.atoms  = Selection.unfold_entities(self.structure[0], "A" )
			# instantanize NeighborSearch
			self.ns = NeighborSearch(self.protein_ca_atoms)
			for residue in self.protein_ca_atoms: 
				residue_attributes = residue.get_full_id()
				if residue_attributes[2] == self.chain :    # check if the residue belongs to the specified chain
					close_residues = self._definepatches(residue, energy_cutoff=energy_cutoff, dist_cutoff=dist_cutoff)
					# Initialize the 'clusters_pre_list'  if it is empty
					if len(clusters_pre_list) == 0 and len(close_residues) != 0 : 
						clusters_pre_list.append(close_residues)
					
					elif len(close_residues) != 0:   # 'close_residues' should not be empty
						# if 'clusters_pre_list' is already initialized, iterate over the list of clusters 
						# to check if any atom of  'close_residues' belongs to one of the clusters
						for cluster in clusters_pre_list : 
							to_add_to_cluster = [ close_residues for aa in close_residues if aa in  cluster ]
							if to_add_to_cluster != []:
								for element in to_add_to_cluster[0]:
									if element not in cluster:
										cluster.append( element )
								break
						else: 
							# when the loop exhausts iterating clusters_pre_list.
							# if neither of the atoms in 'close_residues' belongs to 
							# one of the clusters, then assign all the atoms of 'close_residues'
							# to a new cluster.
							clusters_pre_list.append(close_residues)
			self.clusters_pre_list = clusters_pre_list
		else: 
			raise TypeError ('You have not provided a Foldx alanine scanning file')

	def _definepatches(self, residue, energy_cutoff=2.0, dist_cutoff = 6):	
		"""
		returns a list of of residues which are close and satisfy the energy 
		criteria for being a hotspot.
		"""
		close_residues = self._proximalResidues(residue, distance_cutoff=dist_cutoff)    # what are the residues
		is_close_residue_a_hotspot_table = []
		for close_residue in close_residues : 
			is_close_residue_a_hotspot = self._evaluateEnergy(close_residue, energy_cutoff=energy_cutoff) 
			if not is_close_residue_a_hotspot == None :
				is_close_residue_a_hotspot_table.append(is_close_residue_a_hotspot)
		return is_close_residue_a_hotspot_table
					 
	def _proximalResidues(self, atom, distance_cutoff = 6): 
		"""
		Returns a list of residues which CA atom is close by at least 
		a 'distance_cutoff' to another amino acid.
		"""
		atom_id = atom.get_full_id()[3][1]
		close_atom_list = []
		close = self.ns.search(atom.coord, distance_cutoff )
		for closatom in close:
			close_atom_id = closatom.get_full_id()[3][1]
			if atom_id != close_atom_id:
				close_atom_list.append(closatom)
		return close_atom_list

	def _evaluateEnergy(self, atom, energy_cutoff=3): 
		""" Check if an atom corresponds to a residue with
		high energy or Alanine. return Boolean value.
		"""
		atom_res_id = str( atom.get_full_id()[3][1] )
		
		try: 
			isInTable = self.dic[int(atom_res_id)]
		except: 
			pass 
		
		if isInTable :
			close_residue_name = isInTable[0]
			close_residue_energy = isInTable[1]
			if  close_residue_energy >= energy_cutoff:
				return  atom_res_id, isInTable
		else: 
			return False

	def formatClusters(self):
		"""
		returns a list of 'resid,cluster,chain' for each amino acid
		and a dictionary 'cluster_dic' containing the cluster ID as 
		keys and a tuple as a value contaning the cluster size, 
		volume, and energy density.
		"""
		for chain_cluster in self.clusters_pre_list: 
			cluster_dic = {}
			residue_in_cluster_dic = {}
			residues_to_cluster = {}
			myClusterList = [elem for elem in self.clusters_pre_list  if len(elem) > 3.0   ] 
			all_rsidue_identifiers = [key for key in self.dic.keys()  ]
			# print( all_rsidue_identifiers )
			tags = list( range(1, len(myClusterList)+1 ) )   # ID number of the cluster in the chain
			for tag, cluster in zip(tags, myClusterList ):
				cluster_volume = 0.0
				cumulated_energy = 0.0
				for residue in cluster: 
					cluster_volume += volume[residue[1][0]]
					cumulated_energy += residue[1][1]
					residue_in_cluster_dic[str( residue[0] ) ] = 'c'+str(tag)+residue[1][2]
				per_volume_energy = round( cumulated_energy/cluster_volume, 6 )
				cluster_name =  'c'+str(tag)+cluster[0][1][2]
				cluster_size =  len(cluster) 
				cluster_dic[ cluster_name ] = (cluster_size, round(cluster_volume, 3), per_volume_energy )
			each_residue_to_a_cluster = {}
			output_residue_to_cluster = []
			for residue_id in all_rsidue_identifiers: 
				try :
					each_residue_to_a_cluster[str(residue_id)] = residue_in_cluster_dic[  str(residue_id) ]
					line =  ','.join( [str(residue_id), each_residue_to_a_cluster[str(residue_id)], self.chain] ) 
					output_residue_to_cluster.append( line )
				except: 
					pass # do nothing if the residue is not in a cluster
			
			# Sanity check 	 
			if output_residue_to_cluster : 
				return output_residue_to_cluster, cluster_dic
			else:
				return None   

	def outputToFiles(self, residue_to_cluster, clusters, suffix='Hotspots', path='./', overwrite=True):
		"""
		'residue_to_cluster' is the first list returned by formatClusters
		and 'clusters' is the dictionary returned by formatClusters
		suffix is used to tag the output
		"""
		output_res_clusters = path+'/'+suffix+'_residueClus.csv' 
		output_clusters = path+'/'+suffix+'_clusters.csv'
		if overwrite == True: 
			operation = 'w'
			outputheader1 = ['res_ID', 'Cluster_ID', 'chain']
			outputheader2 = ['Cluster_ID', 'size', 'volume', 'density(kcal/mol/A**3)']
			with open( output_res_clusters , operation) as file1: 
				file1.writelines( ','.join( outputheader1  )+'\n' )
				for line in residue_to_cluster: 
					file1.writelines( str(line)+'\n' )

			with open(output_clusters, operation) as file2: 
				file2.writelines( ','.join( outputheader2  )+'\n' )
				for item in clusters: 
					file2.writelines( ','.join( [item, str(clusters[item][0]) , str(clusters[item][1]), str(clusters[item][2]) ] )+'\n' )

		else:
			operation = 'a'
			with open(output_res_clusters, operation) as file1: 
				for line in residue_to_cluster: 
					file1.writelines( str(line)+'\n' )

			with open(output_clusters, operation) as file2: 
				for item in clusters: 
					file2.writelines( ','.join( [item, str(clusters[item][0]) , str(clusters[item][1]), str(clusters[item][2]) ] )+'\n' )	



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=" A tool to detect the 3D hotspot patches of a protein.")
	# add long and short argument
	parser.add_argument("--pdb", help="PDB file")
	parser.add_argument("--ALAscan", help="Foldx Alanine scan file")
	parser.add_argument("--suffix", help="A tag for the output files. Default = 'Hotspots'.")
	parser.add_argument("--DistanceCutoff", help="Cuoff distance to define a contact between CA atoms. Default = 6 Angstroms")
	parser.add_argument("--EnergyCutoff", help="Energy cutoff above which a residues is considered as a hotspot. Default = 2.0 kcal/mol  ")
	parser.add_argument("--output", help="Path to the output folder")
	args = parser.parse_args()

	# Check if the required arguments are filled
	assert args.pdb != None, 'You must provide a PDB structure'

	assert args.ALAscan != None, 'You must provide an Alanine scanning Foldx file.'

	if args.DistanceCutoff != None:
		distance = float(args.DistanceCutoff)
	else: 
		distance = 6.0 

	if args.EnergyCutoff != None: 
		energy = float(args.EnergyCutoff)
	else: 
		energy = 2.0

	if args.output != None: 
		output = args.output
	else: 
		output = '.'

	if args.suffix != None: 
		suffix = args.suffix
	else: 
		suffix = 'Hotspots'

	# output Verbosity 
	output_file1  = output+'/'+suffix+'_residueClus.csv'
	output_file2  = output+'/'+suffix+'_clusters.csv'

	print("""
	Calculating 3D hotspots: 
				PDB:                         {0}
				Alanine scanning file:       {1}
				CA distance cutoff:          {2}
				Energy cutoff:               {3}
				Output folder:               {4}
				Residues to clusters file:   {5}
				Hotspots output:             {6}

	""".format(args.pdb, args.ALAscan, distance, energy, args.output, output_file1, output_file2 ) )

	# Workflow 
	myala = alaSCanAnalysis(args.pdb)
	mydics = myala.allChainsClusters(args.ALAscan)

	pymoldic = {}
	for chain in mydics: 
		index_of_chain = list(mydics.keys()).index(chain)
		if index_of_chain == 0: 
			overwrite = True
		else: 
			overwrite = False
		myala.DefinePatches(chain, mydics[chain], energy_cutoff=energy,  dist_cutoff=distance) 
		returned_result =  myala.formatClusters() 
		try: 
			myala.outputToFiles(returned_result[0], returned_result[1], suffix=suffix, path=output, overwrite= overwrite)
		except: 
			print("No clusters are found for chain {}".format(chain))


