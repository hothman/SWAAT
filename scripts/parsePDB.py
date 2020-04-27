#!/usr/bin/python3
__author__ = "Houcemeddine Othman"

# Check for that biopython is installed
try:
    from Bio.PDB.PDBParser import PDBParser 
except ImportError:
	raise ImportError('Biopython is a required package')
import warnings   
import sys
import re
import subprocess
from shutil import which
from Bio import pairwise2
import getopt
import csv
import argparse

# check if you have python 3
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

# initiate the parser
parser = argparse.ArgumentParser(description=" A script to map the residue \
IDs in pdb file to the protein sequence coordinates. \
\nCalculates per amino acid SASA, secondary structure and hydrogen bond number ")
# add long and short argument
parser.add_argument("--pdb", help="Input pdb file")
parser.add_argument("--fasta", help="Reference sequences in FASTA format ")
parser.add_argument("--output", help="Path to folder of the output files ")

# read arguments from the command line
args = parser.parse_args()
if not args.pdb:  
    print("A pdb file is required \n Use --help for more details")
    sys.exit(1)
if not args.fasta:  
    print("A fasta file is required \n Use --help for more details")
    sys.exit(1)
if not args.output:  
    output = './'
else:
	output = args.output

# three letters to one letter code of amino acids
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
"""
Reference SASA for extemded conformations (ALA-X-ALA) of amino acids, calculated
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
		   	
class ParsePDB:
	"""handle the PDB files and fill and update the properties"""
	def __init__(self, pdb_file):
		self.pdb_file = pdb_file
		parser = PDBParser( QUIET=True )
		self.structure = parser.get_structure('S',self.pdb_file )
		self.properties = self.getChainSeqResids()

	def getChainSeqResids(self):	
		"""
		Get chain ID, sequence and residue ID of each monomer in the 
		stucture in a list of dictionaries
		"""
		my_structure = self.structure
		if len(my_structure) >  1 :
			raise ValueError('you have multiple models in the PDB file !')
		for model in self.structure :
			monomers = []
			for chain in model: 
				chain_id =  chain.id 
				sequence=''
				resids = []
				for residue in chain: 
					resname = residue.get_resname()
					AA = amino_acids[resname]
					sequence=sequence+AA
					resids.append( residue.id[1] )
				monomer = { 'chain':chain_id, 'seq':sequence, 'ids':resids, 'posinref':[] }
				monomers.append(monomer)
		return monomers

	def sasa(self, chain):
		""" Calculates the fraction of exposed surfqace for each 
		amino acid in the pdb """
		if which("freesasa") == None:
			sys.exit("freesasa is not in the path")
		else:
			constructed_cmd = "freesasa --format seq -n 100 --select "+ " ' anything, chain "+chain+ "' " +self.pdb_file
			output = str(subprocess.check_output(constructed_cmd, shell=True) ) 			
			SASA = re.findall("\d+\.\d+", output)
			letter_sequence = (re.findall("ALA|ASN|ARG|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|TRP", output) )
			ratio_values = []
			sasa_values = []
			for sasa, amino_acid  in zip( SASA, letter_sequence ):
				ratio = float(sasa)/SASAref[amino_acid]
				ratio_values.append(round(ratio,2) )
				sasa_values.append( float(sasa) )
		return sasa_values, ratio_values

	def isInSeq(self,  seq_table  ):
	    """
	    search for returns to monomers and updates the properties list
	    Works on monomers and multimers  """
	    for seq in seq_table:
		    isProt1 = re.match('^[AERTYIPQSDFGHKLMWCVNaertyipqsdfghklmwcvn]+$', seq['sequence'])
		    if bool(isProt1 )  == True :
		    	for chain in self.properties: 
		    		isProt2 = re.match('^[AERTYIPQSDFGHKLMWCVNaertyipqsdfghklmwcvn]+$', chain['seq'] )
		    		if bool(isProt2 )  == True : 
		    			is_sub_string = seq['sequence'].find(chain['seq'])
		    			if is_sub_string != -1:
		    				length_reference = len(chain['seq'])
		    				length_chain = len( seq['sequence'] )		
		    				start = is_sub_string
		    				end =  length_chain+9 
		    				index_structure = list( range(start, end ) )
		    				chain['posinref'] = index_structure
		    				chain['header'] = seq['header']
		    else:
		    	sys.exit("One of the sequences contains non standard amino acids")
	    for chain in self.properties:
	    	if chain['posinref'] == []:
	    		raise MissingDataError( str(chain['chain']), " Does not match any elements in the sequence table.\n \
	    		Check for mismatches between the reference and the PDB sequences or for missing segments in the PDB file")
	    		sys.exit()

	def getSsHbonds(self, chain):
		""" Secondary structure and hydrogen bond computation parser for stride output"""
		if which("stride") == None:
			sys.exit("stride is not in the path")
		else: 
			# get the secondary structure
			constructed_cmd = "stride -f "+self.pdb_file+" -r"+chain+" -h |grep ASG"
			output = str(subprocess.check_output(constructed_cmd, shell=True) )
			byline_split = output.split('\\n')
			sec_structure = []
			for line in byline_split:
				byfield_split =  re.compile("\s+").split(line) 
				try : 
					sec_structure.append(byfield_split[5])
				except:
					pass
			# get H-bonds
			constructed_cmd = "stride -f "+self.pdb_file+" -r"+chain+" -h |grep 'ACC\|DNR'"
			output = str(subprocess.check_output(constructed_cmd, shell=True) )
			byline_split = output.split('\\n')
			H_bonds=[]
			for line in byline_split:
				byfield_split =  re.compile("\s+").split(line) 
				try : 
					H_bonds.append( int(byfield_split[3]) )
				except:
					pass
		# count the number of H bonds per aa
		hbond_count = []
		for monomer in self.properties :
			if monomer['chain'] == chain:
				for res_id in monomer['ids'] : 
					count = H_bonds.count(res_id)
					hbond_count.append( count )
		return hbond_count, sec_structure

	def fillProperties(self):
		"""  Add the missing data to self.properties """
		for monomer in self.properties: 
			sasa_chain, ratio_chain = self.sasa( monomer['chain'] ) 
			monomer['sasa'] = sasa_chain
			monomer['sasa_ratio'] = ratio_chain
			hbonds, s_structure = self.getSsHbonds(monomer['chain'])
			monomer['SS'] = s_structure
			monomer['hbonds'] = hbonds

	def propertiesToTsv(self, output):
		"""  Output self.properties to tsv file """
		for chain in self.properties: 
			outputfile = output+'/'+chain['header'][1]+'_'+chain['chain']+'.tsv'
			basename = self.pdb_file.split('/')[-1]
			with open( outputfile , 'w') as file:
				file.writelines( '#'+' '.join(chain['header'] )+' '+basename+'\n' )
				TSVhead = 'AA\tID\tIDref\tsasa\tsasa_ratio\tSS\tn_HBonds\n'
				file.writelines( TSVhead )
				for seq, ids, posinref, sasa, sasa_ratio, ss, h_bonds in zip( chain['seq'], chain['ids'], chain['posinref'], chain['sasa'], chain['sasa_ratio'], chain['SS'], chain['hbonds'] ) :
					file.writelines( seq+'\t'+str(ids)+'\t'+str(posinref)+'\t'+str(sasa)+'\t'+str(sasa_ratio)+'\t'+ss+'\t'+str(h_bonds)+'\n' )
			file.close()

class ParseFASTA:
    """ Exxtract the header and the sequence
    Assumes that the header is arranged as follows
    UNIPROTID|GeneName|ENSEMBLgeneID|TranscriptID
    """
    def __init__(self, filename):
        self.filename=filename

    def readFASTA(self):
        self.fastafile = open(self.filename, 'r')
        self.seqlist = self.fastafile.read().split('>')[1:]

        #if len(self.seqlist) > 1 :
        #    sys.exit("The Fasta file contains more than one sequence")
        my_seqs = []
        for mySeq in self.seqlist:
        	li = mySeq.split('\n', 1)
        	header = li[0].split("|")
        	seq =  re.sub('\W+','', li[1] )
        	seq.strip()
        	my_seqs.append( { 'header': header, 'sequence':seq.strip()})
        self.my_seqs = my_seqs
        return self.my_seqs

def alignProteinseq( seq1, seq2):
	matrix = matlist.blosum62
	isProt1 = re.match('^[AERTYIPQSDFGHKLMWCVNaertyipqsdfghklmwcvn]+$', seq1)
	isProt2 = re.match('^[AERTYIPQSDFGHKLMWCVNaertyipqsdfghklmwcvn]+$', seq2)

	if bool(isProt1 )  == True and bool(isProt2 )  == True :
		alignments = pairwise2.align.globalxx(seq1, seq2)
		# print(alignments)
		return alignments
	else:
		sys.exit("One or both sequences contain non standard amino acids")

##########################################
			#The workflow

# Parse the reference fasta file
my_sequence = ParseFASTA(args.fasta)
fileSeqs =  my_sequence.readFASTA()
# Parse the protein PDB file
my_pdb = ParsePDB(args.pdb)
# generate themapping betwenn PDB and fasta
my_pdb.isInSeq( fileSeqs )
# Calculate the per AA sasa and the hbonds and SS
my_pdb.fillProperties()
# outpu to tsv files
my_pdb.propertiesToTsv(output)