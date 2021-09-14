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

## Running the workflow 

Options of the workflow can be listed using the following command.

```bash
nextflow run prepare_db.nf --help 
```

The user is required to provide the list of Uniprot IDs corresponding to the genes to annotate as well as their molecular structures in PDB file format.  

* The structure must have a single chain.
* A structure must contain a chain ID. 
* Heteroatoms are not taken into account
