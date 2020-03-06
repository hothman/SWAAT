# usage: 
# Rscript fetch_seq.R UNIPROTLIST OUTPUTDIR
# UNIPROTLIST is text file containing a header Identifier in the first line
# and one Uniprot ID per line
#

library(biomaRt)
library(tidyverse)

args <- commandArgs(TRUE)
inputList <- as.character(args[1])
outFolder <- as.character(args[2])
setwd(getwd() )

# get the sequences from refseq IDs
getLonguest <- function(x) {
  refseq<-x[3]
  protID<-x[1]
  myseq<- getSequence(id=refseq, mart = ensembl, type='refseq_mrna', seqType = 'peptide')
  return( as.character(myseq$peptide)  )
}

# select the loguest isoform and format to FASTA
output2Fasta <- function(protID,DF, outputfolder) {
  # ProtID is a uniprot identifier, DF: a data frame containing isoforms of many proteins
  allIsoforms <- filter(DF, DF$uniprotswissprot== protID)
  maxlength <- filter(allIsoforms, nchar(allIsoforms$Sequences ) == max(nchar(allIsoforms$Sequences)) ) 
  selected <- maxlength[1,] 
  seq <- str_remove(as.character( selected$Sequences ), '[*\t\r\f]')
  Uniprot <- selected$uniprotswissprot
  hgnc <- selected$hgnc_symbol
  ensembl <- selected$ensembl_gene_id 
  refseq_mrna <- selected$refseq_mrna
  fileoutput <- paste(outputfolder,"/",hgnc,".fa", sep="")
  #format to fasta string
  fasta <- paste(">",Uniprot,"|",hgnc,"|",ensembl, "|", refseq_mrna , '\n', seq, sep="")
  write( fasta, fileoutput )
  return(fasta)
}

# Use grch37 genome assembly
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                  path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
print(listDatasets(ensembl)[26,2])

# read the uniprot ID
protID <- as.vector( read.csv( inputList ))

# Fetch the sequences
BM<- getBM(filters='uniprotswissprot', 
      attributes = c('uniprotswissprot', 'ensembl_gene_id', 'refseq_mrna', 'hgnc_symbol'), 
       mart = ensembl, values=protID)

# discard entries without refseqID
BM <- filter(BM, BM$refseq_mrna != "")
UniqID <- unique(BM$uniprotswissprot)

Sequences <- apply(BM, 1, getLonguest)
BM$Sequences <- Sequences

dir.create(outFolder )

myFasta <- lapply(UniqID  ,output2Fasta , DF=BM, outputfolder=outFolder )