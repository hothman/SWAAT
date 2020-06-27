#!/usr/bin/python3
""" This script takes the variant file relative to each gene (tsv), formats
        and generates the input for in silico mutagenesis
usage:
./parseVar.py -i variants.dat -o outputfolder -m genemap.tsv -f gene_data.tsv
"""
__author__      = "Houcemeddine Othman"


import sys
import csv
import os
from shutil import rmtree
import re
import getopt

import getopt

# check if you have python 3
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

# parse the arguments list in getopt() style
"""
    -i input variants file
    -o output folder
    -f path to chain info file
    -m mapping file
"""
args, out=getopt.getopt(sys.argv[1:], 'i:f:o:m:' )
argdic={}

for elem in args:
    argdic[elem[0]] = elem[1]

def checkfiles(geneData, inputvar, mapfile):
    '''
    Checks if the input files are accessible
    '''
    if os.path.isfile(geneData) == False :
        print( 'Gene data file is not found')
        sys.exit(1)
    elif os.path.isfile(inputvar) == False :
        print( 'Variants list file is not found')
        sys.exit(1)
    elif os.path.isfile(mapfile) == False  :
        print( 'Mapping file is missing')
        sys.exit(1)

# Check that the input list matches
def readfile(inputvar, geneData):
    """
    Check that the input list matches the
    database gene names (exit 1 if not)
    """
    dic={}
    for line in inputvar:
        theLine = line.split("\t")
        if theLine[0] not in geneData:
            print(theLine[0], "does not seem to be an ADME gene")
            sys.exit(1)
        elif theLine[0] in dic.keys():
            dic[ theLine[0] ] =  dic[ theLine[0] ]+" "+theLine[1].replace('\n','')
        else:
            dic[ theLine[0] ] = theLine[1].replace('\n','')
    print('***** All variants match the ADME gene list *****')
    return dic

def formatDB( geneData):
    """
    Converts the database file to dict and pull
    the gene list
    """
    with open( geneData ) as tsvfile:
      reader = csv.DictReader(tsvfile, dialect='excel-tab')
      mygenes= list( reader )
      #[ geneList.append(row['ID']) for row in reader ]
    geneList =[]
    [ geneList.append(row['ID']) for row in mygenes ]
    geneinfo = {}
    for elem in mygenes:
        geneinfo[ elem['ID'] ]=elem['chain'].split(',')
    return geneinfo, geneList

def createfolders(outputFolder, inputvar):
    '''
    Reads the variants dictionary outputted from
    'readfile' function and generates folders
    '''
    # creating the parent folder
    try:
        os.mkdir(outputFolder)
    except FileExistsError:
        rmtree(outputFolder)
        os.mkdir(outputFolder)

    #creating the child folder
    for elem in inputvar:
        # print( elem['ID'] )
        try:
             os.mkdir( outputFolder+'/'+elem )
        except FileExistsError:
             os.rmdir( outputFolder+'/'+elem )
             os.mkdir( outputFolder+'/'+elem )
        allVariants = inputvar[elem].split(' ')

def parseMap(pathmapfile):
    '''
    Takes a mapping tsv file havin three fields. The
    name of the genes preceeds the table and the entries
    are concatenated
    output: a dictionary (genemap_dic) containing the gene ID as key
    and a list of the mapping data as values
    example:
    {'gene2': [(1, 22, 'G'), (2, 23, 'V'), (3, 24, 'E')],
    'gene3': [(1, 22, 'G'), (2, 23, 'V'), (3, 24, 'E')],
    'gene1': [(1, 22, 'G'), (2, 23, 'V'), (3, 24, 'E')]}
    '''
    with open(pathmapfile, 'r') as infile:
        data = infile.readlines()
    # each block corresponding to a gene
    # is proceeded as
    dic={}
    index=[]
    for i ,elem in enumerate(data):
        if '#' in elem:
            index.append(i)
    chunks=[]
    for elem in reversed(index):
        chunks.append( data[elem:])
        del data[elem:]
    # process data
    genemap_dic = {}
    for gene in chunks:
         residues = {}
         key=gene[0].replace('#','').strip()
         for AA in  gene[1:]:
             if len( AA.split('\t') ) != 3:
                 pass
             else :
                 splitted = AA.split('\t')
                 # position0: aa id in pdb model,
                 # position1: aa id in reference sequence
                 # position2: amino acid in ref seq
                 residues[ int(splitted[0]) ] = ( int(splitted[2]), splitted[1] )
         genemap_dic[key] = residues
    return(genemap_dic)

def createIndiv(variants, geneinfo, genemap_dic):
    '''
    generates a dictionary containing formatted strings, ready to be
    written as individual filesself.
        variants: output of readfile() function
        geneinfo: output dictionary of formatDB()
        genemap_dic:
    '''
    output = {}
    # each elem is a gene
    for elem in variants:
        elem_map = genemap_dic[elem]       # read the map for the gene elem
        elem_chain = geneinfo[elem]        # read the chains ID for the gene elem

        # parsing each haplotype alone
        construction_list = []
        for var in variants[elem].split(' '):
            coor_ref = [ int(number) for number in re.findall(r'\d+', var)  ]
            aa_ref = re.findall(r'[XAERTYIPQSDFGHKLMWCVNxaertyipqsdfghklmwcvn]', var)

            # construct the formatted string
            table4haplotypes = []
            for id, aa in zip(coor_ref, aa_ref):
                for chain in elem_chain:
                    construction =  elem_map[id ][1]+chain+str(elem_map[id ][0])+aa
                    table4haplotypes.append( construction )
            construction_list.append(','.join(table4haplotypes)+';')
        output[elem] = construction_list
    return output

def writeOutput(output, outputFolder):
    '''
    Writes the formatted string to the outputFolder
    '''
    for gene in output:
        i = 1
        for haplotype in output[gene]:
            filename = str(gene)+'_h'+str(i)+'.txt'
            #subfoldername = 'haplotype'+str(i)
            PATH2OUTPUT='/'.join([outputFolder, gene, filename])
            PATH2FOLEROUTPUT = '/'.join([outputFolder, gene])
            # try:
            #     os.mkdir(PATH2FOLEROUTPUT)
            # except:
            #     os.rmdir(PATH2FOLEROUTPUT)
            #     os.mkdir(PATH2FOLEROUTPUT)
            with open(PATH2OUTPUT, 'w') as text_output:
                text_output.write(haplotype+'\n')
            print("Writing %s to %s " % (haplotype, PATH2OUTPUT ))
            path = os.getcwd()
            i += 1

def runPipeline(input, chain_info ,map_file, outputfolder ):
    # check input file
    checkfiles(chain_info, input, map_file )
    # read the file containing the protein chain data
    chains, genes = formatDB( chain_info )
    # read input file containing variants
    file_name = open ( str( input ), 'r' )
    variants = readfile( file_name, genes)
    # pqarse the map file
    parsed_data = parseMap( map_file )
    #create the formatted strings dictionary
    formattedstr = createIndiv(variants, chains, parsed_data )
    # create containing createfolders
    createfolders(outputfolder, variants)
    # write to folders
    writeOutput( formattedstr, outputfolder )
    print('***** Finished writing to %s ***** ' % outputfolder)

#########################################################################
# generate the input files

chain_info = argdic['-f']
map_file = argdic['-m']
inputfile = argdic['-i']
output = argdic['-o']
runPipeline(inputfile, chain_info, map_file, output )
