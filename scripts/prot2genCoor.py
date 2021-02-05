#!/usr/bin/python3
"""
Usage:
prot2genCoor -i<input-fasta> -d <headers> -o <output-file>

The script does not tolerate multiple FASTA files
"""
__author__ = "Houcemeddine Othman"

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
args, out=getopt.getopt(sys.argv[1:], 'i:o:f:' )
argdic={}
for elem in args:
    argdic[elem[0]] = elem[1]

# Check the input fasta
try:
    fastaInput = argdic['-i']
except:
    sys.exit('No FASTA file is provided. We stop here.\
    \nUsage:\nprot2genCoor -i<input-fasta> -o <output-file> -f <output_fasta> ')

# Check the output
try:
    output = argdic['-o']
except:
    sys.exit('No output is provided. We stop here. \
    \nUsage:\nprot2genCoor -i<input-fasta> -o <output-file> -f <output_fasta> ')

try:
    fasta_output = argdic['-f']
except:
    sys.exit('No output is provided for the refseq protein file. We stop here. \
    \nUsage:\nprot2genCoor -i<input-fasta> -o <output-file> -f <output_fasta> ')


## Process FASTA file
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

        if len(self.seqlist) > 1 :
            sys.exit("The Fasta file contains more than one sequence")

        mySeq = self.seqlist[0]
        #print(mySeq)
        li = mySeq.split('\n', 1)
        header = li[0].split("|")
        seq =  re.sub('\W+','', li[1] )
        self.header = header
        self.seq = seq
        return self.header, self.seq
 
class Transvar():
       """Run transvar for each amino acid position and 
       generate a tsv table containing 
       chromosomeID, HGCN name, amino acid, amino acid number
       starting and ending positions of the codon and 
       both cDNA and gDNA codents"""
       def __init__(self, header, sequence, outputfile, fasta_output):
           self.header = header 
           self.sequence = sequence
           self.outputfile = outputfile
           self.fasta_output = fasta_output
           print(self.header)
       
       def getRegexp(self):
           """ process the attributes and generate the output """

           coord_list = []
           transcript_name = self.header[3]
           print("Running transvar for transcript "+transcript_name+" of "+self.header[1]+" ...")
           geneName= self.header[1]

           aa_index = 1
           output = "initial output"

           ref_seq_aa_sequence = ""
           while output != "":
            try:
              constructed_cmd = "transvar  panno -i '"+str( geneName+":p."+str(aa_index) )+"' --refseq --refversion hg19|grep "+transcript_name
              output = str(subprocess.check_output(constructed_cmd, shell=True) )
              regexp1 = re.compile(r'chr\d+:g.\d+_\d+/c.\d+_\d+/p.\d+\w') 
              gen_coor = regexp1.findall(output)
              regexp2 = re.compile(r'gDNA_sequence=\w{3}')
              gDNA = regexp2.findall(output)
              regexp3 = re.compile(r'cDNA_sequence=\w{3}')
              cDNA = regexp3.findall(output)
              regexp4 = re.compile(r'protein_coding[)]\\t\w+\\t')
              gene_hgcn = regexp4.findall(output)
              amino_acid = gen_coor[-1][-1] 
              if amino_acid in "QWERTYIPASDFGHKLXCVNM":
                ref_seq_aa_sequence = ref_seq_aa_sequence + gen_coor[-1][-1]
              reg1 = re.compile(r'\d+')
              exp1 = reg1.findall(gen_coor[0])
              start_genomic_position = exp1[1]
              end_genomic_position = str (int( start_genomic_position ) + 2 )
              AA_position = exp1[-1]
              chromosome = exp1[0]
              AA_ref = gen_coor[0][-1]
              gene_hgcn = gene_hgcn[0].split('\\t')[1]
              gDNA_codon = gDNA[0].replace("gDNA_sequence=", "")
              cDNA_codon = cDNA[0].replace("cDNA_sequence=", "")
              coord_list.append([chromosome, gene_hgcn, AA_ref, AA_position, start_genomic_position, \
                 end_genomic_position, gDNA_codon, cDNA_codon] )

            except: 
              break
            aa_index = aa_index +1 

           with open( self.outputfile , 'w') as file:
            file.writelines("chr\tgene_hgcn\tAmino_acid\tAA_position\tgen_start\tgen_end\tgDNA_codon\tcDNA_codon\n")
            file.writelines('\t'.join(i) + '\n' for i in coord_list) 

           with open( self.fasta_output , 'w') as refseq_fasta: # output the refseq protein sequence to a fasta file 
            refseq_fasta.writelines(">"+self.header[0]+"|"+ self.header[1]+"|"+ self.header[2]+"|"+self.header[3]+"\n")
            refseq_fasta.writelines(ref_seq_aa_sequence)

           print("Finished successfully")   

class CheckConsitency(object):
   """This class is used to check the consistencey of the sequence data between the Uniprot and refseq"""
   def __init__(self):
     pass
       
####################
# Run the workflow #
#################### 
if __name__ == "__main__":
     Fastafile= ParseFASTA(fastaInput)                         # createarseFASTA() object
     header, sequence =Fastafile.readFASTA()                   # extract header and sequence
     obj = Transvar(header, sequence, output, fasta_output)    # create transvar object
     obj.getRegexp()                                           # run transvar and generate tsv table
