# SWAAT

![Drag Racing](logo_SWAAT.png)

## Preparing dependencies for the auxiliary workflow 

SWAAT is composed of the main workflow that annotates a list of 36 ADME genes and an auxilliary workflow that could be used to prepare the database for any gene other than from the list of 36 ADME members.  

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



