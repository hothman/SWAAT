# Auxiliary workflow to run SWAAT annotation on non-ADME genes

The auxiliary workflow will help in building the database required to annotate other non ADME genes that are not available within the default settings that we have provided with SWAAT. 

## Dependencies 

* Python 3.7.4
* [Nextflow](https://www.nextflow.io/) 20.10.0
* [Transvar](https://transvar.readthedocs.io/en/latest/index.html) 2.5.9
* Biopython  1.76
* [PRODRES](https://github.com/ElofssonLab/PRODRES)
* [FoldX](http://foldxsuite.crg.eu/)
* [ENCoM](https://github.com/NRGlab/ENCoM)
* [freesasa](https://freesasa.github.io/) 2.0.3

All the dependencies are open source. The user needs only to obtain a free licence for [FoldX](http://foldxsuite.crg.eu/). We recommend installing the dependencies using conda virtual environment.  We also recommend using the [Mamba](https://github.com/mamba-org/mamba) package manager which will speed up the building process. 

A special consideration must be paid to the installation of ENCoM. ENCoM can be obtained from its [git repository](https://github.com/NRGlab/ENCoM). The `build_encom` code systematically returns an exit status of 1 which creates an issue when running within  a Nextflow process. To overcome the problem, you need to modify line 370 from `return(1);` to `return(0);` in `src/build_encom.c`. Thereafter you can compile the set of codes as following: 

```sh
cd ENCoM/build/ \
	&& make -f Makefile.Unix
```



## Input data

The user is required to provide:

*  A list of Uniprot IDs corresponding to the genes to annotate
* Molecular structures in PDB file format corresponding to the proteins encoded by the genes to annotate. 

### Providing the list of genes

The Uniprot list must be arranged in a text file **containing no header** as one accession per line. No header or hanging empty lines should be provided in the file. A typical layout is the following. 

```
O75795
P20813
P33261
Q16348
```

### Preparing the PDB structures

Follow these guidance to prepare the structures. 

* The structure must have a single chain.
* The PDB file must be named as follows `<UNIPROTID>.pdb` (e.g. P04637.pdb )
* A structure must contain a chain ID. 
* Heteroatoms are not taken into account
* No discontinuities are allowed in the protein chain. 

The [MMTSB tool set](https://github.com/mmtsb/toolset) is could be very handy to prepare the protein structure.

### Calculation of the PSSMs

The process that runs this task uses the standalone version of [PRODRES](https://github.com/ElofssonLab/PRODRES). The user must set the corresponding options in the command line or in the workflow file. However, if PRODRES is not installed, the generation of the PSSMs will not take place, However, Nextflow won't complai about this. Notice also that the PSSMs are mandatory when using the main workflow of annotation. We have made this choice because the installation of PRODRES may be non trivial for the user. SWAAT and the auxiliary workflow **use the protein sequences of the Refseq database**. 

However, the user can download the FASTA Refseq sequences of the proteins, modify their header lines and and submit them to the [server version of PRODRES](https://prodres.bioinfo.se/).  However, we recommand running the auxiliary workflow first which will prepare the FASTA sequences for you, and then submit them to the server as a multi sequence FASTA file/string.  

### FTMap input data 

FTMap is a lool ([doi.org/10.1038/nprot.2015.043](doi.org/10.1038/nprot.2015.043) ) used in the Identification of binding hot spots. This is the only step that will require a manual processing by the user. However, if the gene is not known to bind drugs or small ligands, it won't be necessary to go through this step. The user needs to submit the structure to the  [STMap server](http://ftmap.bu.edu/login.php) and then download the 'hbonded' and the 'nonbonded' raw output output files to a clean location with the **gene symbol as the name for the directory**, The files should be tagged "nonbonded" and "hbonded" strings (e.g. `mage6_hbonded.rawextract` and `mage6_nonbonded.rawextract`) The location of the parent directory, i.e. containing the gene symbol named folders, can then be specified in the `ftmap` option.

## Running the workflow 

Options of the workflow can be listed using the following command.

```bash
nextflow run prepare_db.nf --help 
```

A typical running of the workflow is as follows: 

```
nextflow run prepare_db.nf --protlist /path/to/Uniprot_list.txt \
	--pdbs /path/to/pdb/files
	--outfolder /path/to/the/outputfolder
```

## Output

The auxiliary workflows generates a series of directories, each containing data that are used by the main workflow to speed-up the annotation process. 

* `sequences`: Contains the FASTA sequences extracted from Refseq and Uniprot. 
* `uniprot2PDBmap` contains the mapping of the reference sequence protein coordinates to the PDB sequence coordinates. 
* `Seq2Chain` links the gene symbol to Uniprot ID, the chain identifier and the PDB file name.
* `prot_annotation` offers information about the functional contribution of amino acids in the reference sequence. 
* `maps` contains the mapping between the gene coordinates and the protein reference sequence. 
* `hotspots` contains the identified hotspot information relative to the PDB coordinates. 
* `ENCom` contains the normal modes calculated for the reference protein structure.
* `PSSMs` contains the pssm files calculated for the protein using the reference sequence.
* `PDBs` contains the 3D structures in PDB file format.
* `ftmap` contains scored amino acid positions for ligand/drug binding putative site. 

## Example: Melanoma-associated antigen 6 (MAGE-6)

The Uniprot accession for MAGE-6 is [P43360](https://www.uniprot.org/uniprot/P43360). The structure of this protein contains one chain. The structure were generated using homology modeling and provided in the directory `input_example/PDB` as `P43360.pdb`. 

Next you need a text file containing no header in the following format

```
P43360
```

The file can also be found in the `input_example` directory (`example.csv`) . 

The FTMap files were calculated by submitting the PDB structure to the servers which results in two output files. They were renamed to contain the tags "hbonded" and "nonbonded". The raw files were put into a directory `MAGEA6` named after the gene symbol (this is a requirement). The `MAGEA6` directory is a child directory of `input_example/ftmap`. We can then use the option `--ftmap` to point to the location of the ftmap raw files.  

One you have prepared the inputs you can call the auxiliary workflow as follows: 

```bash
nextflow run prepare_db.nf --protlist /home/hothman/SWAAT/auxiliary_wf/input_example/example.csv \
	--pdbs /home/hothman/SWAAT/auxiliary_wf/input_example/PDB \
	--outfolder /home/hothman/new_genes_database
	--ftmap /home/houcem/SWAAT/auxiliary_wf/input_example/ftmap
```

A successful run will create a new folder containing the database directory (`new_genes_database`) that are necessary to annotate the MAGE-6 gene using structural data. The output of the command should give something similar to the following:

```
Number of genes    : 1
Path to PDB files  : /home/houcem/SWAAT/auxiliary_wf/input_example/PDB
Outputdir          : database_new_gene/
List of accessions : P43360
executor >  local (8)
[d9/9641a8] process > GetProteinAnnotationFetchFasta (1) [100%] 1 of 1 ✔
[61/310184] process > GetCoordinates (1)                 [100%] 1 of 1 ✔
executor >  local (9)
[d9/9641a8] process > GetProteinAnnotationFetchFasta (1) [100%] 1 of 1 ✔
[61/310184] process > GetCoordinates (1)                 [100%] 1 of 1 ✔
[42/fb35d6] process > geneToChainMapping (1)             [100%] 1 of 1 ✔
[58/0178c7] process > uniprot2PDB (1)                    [100%] 1 of 1 ✔
[ab/c18a85] process > Hospotislands (1)                  [100%] 1 of 1 ✔
[6a/76b297] process > encomWT (1)                        [100%] 1 of 1 ✔
[04/7ba5b9] process > parseFTMAP (1)                     [100%] 1 of 1 ✔
[a3/4dca9d] process > generate_matrixes                  [100%] 1 of 1 ✔
[d8/c557ae] process > outputPDB (1)                      [100%] 1 of 1 ✔
Completed at: 06-Nov-2021 14:15:35
Duration    : 1m 53s
CPU hours   : (a few seconds)
Succeeded   : 9

```



