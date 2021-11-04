![Drag Racing](logo_SWAAT.png)

# Running SWAAT

First, you need to clone the SWAAT repository from GitHub: 

```sh
git clone https://github.com/hothman/SWAAT.git
```

You need than to change your path to the repository where you can run the workflow from. The main workflow, that runs the annotation process is implemented in `main.nf` file. To list the options and the arguments type the follwowing command in the source directory after  :

```
nextflow main.nf --help

Arguments:
  --dbhome [folder]               Path to database containing the dependency files for annotating the variants (Default False)
  --vcfhome [folder]              Path to folder containing VCF files split by annotated gene (e.g. CYP2D6.vcf) (Default False)
  --outfolder [str]               Where to output the plain text and the HTML report (Default: false)
  --genelist [file]               User can limit the annotation to the list of genes contained in a this text file (one line per gene) (Default False)

Other
  --foldxexe [str]                Specifies the name of the executable of FoldX software (Default foldx)
  --encomexe [str]                Specifies the name of the executable of build_encom (Default build_encom)
  --freesasaexe [str]             Specifies the name of the executable of freesasa software (Default freesasa)
  --strideexe [str]               Specifies the name of the executable of stride software (Default freesasa)
  --rotabase [abs path]           Path to rotabase.txt (Default PATH to foldx)
```

The VCF files to annotate should be specific for each ADME gene in an uncompressed format. The file name for the VCF file must contain the gene symbol in upper case letters (e.g CYP1A2.vcf and TPMT.vcf). A full list of the annotated ADME genes is contained in `./database/gene_list/all_adme_genes.txt`. The path to folder that contains the VCF files is specified in the `--vcfhome` argument. 

A typical run of the workflow could be the following: 

```sh
nextflow run main.nf --dbhome /home/hothman/SWAAT/database/ \
		--vcfhome /home/hothman/SWAAT/vcfs \
		--outfolder ./swaat_out \
		--genelist ./inputexample/gene_list.txt
```

A maintained database can be downloaded with the repository in `./database` folder. The argument `--dbhome` allows to specify the path to the database.



## Preparing dependencies for the auxiliary workflow 

SWAAT is composed of the main workflow that annotates a list of 36 ADME genes and an auxiliary workflow that could be used to prepare the database for any gene other than from the list of 36 ADME members.  

### Dependencies for the main workflow

* Python 3.7.4
* [Nextflow](https://www.nextflow.io/) 20.10.0
* Biopython  1.76
* numpy 1.18.1
* pandas 1.0.4 
* scikit-learn 0.22.2
* scipy 1.4.1 
* [FoldX](http://foldxsuite.crg.eu/)
* [ENCoM](https://github.com/NRGlab/ENCoM)
* [freesasa](https://freesasa.github.io/) 2.0.3
* [Stride](http://webclu.bio.wzw.tum.de/stride/) 
* [DSSP ](https://github.com/cmbi/dssp)  >= 2.2.6

### Dependencies for the auxiliary workflow

* Python 3.7.4
* [Nextflow](https://www.nextflow.io/) 20.10.0
* [Transvar](https://transvar.readthedocs.io/en/latest/index.html) 2.5.9
* Biopython  1.76
* [PRODRES](https://github.com/ElofssonLab/PRODRES)
* [FoldX](http://foldxsuite.crg.eu/)
* [ENCoM](https://github.com/NRGlab/ENCoM)
* [freesasa](https://freesasa.github.io/) 2.0.3

### Instruction for installing the dependencies

All the dependencies are open source. The user needs only to obtain a free licence for [FoldX](http://foldxsuite.crg.eu/). We recommend installing the dependencies for SWAAT using conda virtual environment.  First let's create a virtual environment called `SWAAT` and install all the python dependencies including biopython, numpy, pandas, scikit-learn and scipy. We provided a `yaml` file in the root directory to automate the process. We also recommend using the [Mamba](https://github.com/mamba-org/mamba) package manager which will speed up the building process. 

First, you would have to install [Anaconda](https://www.anaconda.com/products/individual) or [miniconda](https://docs.conda.io/en/latest/miniconda.html) to be able to create and manage virtual environments. 

To create and build an environment called SWAAT use the following command: 

```sh
conda env create -f SWAAT.yml
```

Once the environment is built you can activate it using the following command:

```sh
conda activate SWAAT
```

 [Transvar](https://transvar.readthedocs.io/en/latest/index.html) is not required to annotate the list of the default ADME genes. However you will need the tool if you want to run the auxiliary workflow. 

```sh
pip install transvar
# set up databases
transvar config --download_anno --refversion hg19
# in case you don't have a reference
transvar config --download_ref --refversion hg19
# test
transvar panno -i 'PIK3CA:p.E545K' --ucsc --ccds
```
You can also configure transvar to use the hg38 build of the reference genome.


### Installing dssp

```
git clone https://github.com/cmbi/dssp.git
cd dssp 
./autogen.sh
./configure
make
```

### Preparation of the PDB structures
A structuer is a chain
A structure must contain a chain identifier
Heteroatoms are not taken into account

### Installing ENCoM

ENCoM can be obtained from its [git repository](https://github.com/NRGlab/ENCoM). The `build_encom` code systematically returns an exit status of 1 which creates an issue when running within  a Nextflow process. To overcome the problem, you need to modify line 370 from `return(1);` to `return(0);` in `src/build_encom.c`. Thereafter you can compile the set of codes as following: 

```shell
cd ENCoM/build/ \
	&& make -f Makefile.Unix
```

 



