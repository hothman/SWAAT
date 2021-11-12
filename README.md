![Drag Racing](logo_SWAAT.png)



SWAAT allows for the prediction of the variant effect based on structural properties. SWAAT annotates a panel of 36 ADME genes including 22 out of the 23 clinically important members identified by the PharmVar consortium. The workflow consists of a set of python codes. The execution of these code is managed within Nextflow to annotate coding variants based on 37 criteria. SWAAT also includes an [auxiliary workflow](./auxiliary_wf/README.md) allowing a versatile use for genes other than the ADME members. The auxiliary workflow can be found in the directory `auxiliary_wf` of the repository.

# 1. Installing dependencies 

## 1.1 Dependencies for the main workflow

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

## 1.2 Instruction for installing the dependencies

All the dependencies are open source. The user needs only to obtain a free licence for [FoldX](http://foldxsuite.crg.eu/). We recommend installing the dependencies for SWAAT using conda virtual environment whenever it's possible.  First let's create a virtual environment called `SWAAT` and install all the python dependencies including biopython, numpy, pandas, scikit-learn and scipy. We provided a YAML file in the root directory to automate the process. We also recommend using the [Mamba](https://github.com/mamba-org/mamba) package manager which will speed up the building process. 

First, install [Anaconda](https://www.anaconda.com/products/individual) or [miniconda](https://docs.conda.io/en/latest/miniconda.html) to be able to create and manage virtual environments. 

To create and build an environment called SWAAT use the following command: 

```sh
conda env create -f SWAAT_dependencies.yml
```

Once the environment is built you can activate it using the following command:

```sh
conda activate SWAAT
```



## 1.3 Special consideration for Installing ENCoM

ENCoM can be obtained from its [git repository](https://github.com/NRGlab/ENCoM). The `build_encom` code systematically returns an exit status of 1 which creates an issue when running within  a Nextflow process. To overcome the problem, you need to modify line 370 from `return(1);` to `return(0);` in `src/build_encom.c`. Thereafter you can compile the set of codes as following: 

```shell
cd ENCoM/build/ \
	&& make -f Makefile.Unix
```

# 2. Running SWAAT

First, you need to clone the SWAAT repository from GitHub: 

```sh
git clone https://github.com/hothman/SWAAT.git
cd SWAAT
```

The main workflow, that runs the annotation process is implemented in `main.nf` file. To list the options and the arguments type the following command in the source directory after  :

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

# 3. SWAAT Output

SWAAT generates a CVS file, a HTML file, and a set of directories referring to each annotated gene. The output are complementary. 

## 3.1 CSV output

The CSV file offers richer content than the HTML file and serves for elaborate analysis and filtering methods. Each column provides a piece of information according to structure and sequence-based features.

| Column                                                                                                                                                                                                                                                                                                                                | information                                                                                                                                                                                   |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `gene_name`                                                                                                                                                                                                                                                                                                                           | Gene symbol name                                                                                                                                                                              |
| `wt_res`                                                                                                                                                                                                                                                                                                                              | Residue type in the amino acid reference sequence                                                                                                                                             |
| `mut_res`                                                                                                                                                                                                                                                                                                                             | Amino acid type encoded by the genetic variant in the provided VCF file                                                                                                                       |
| `position`                                                                                                                                                                                                                                                                                                                            | Position in the amino acid reference sequence                                                                                                                                                 |
| `chain`                                                                                                                                                                                                                                                                                                                               | Chain identifier in the PDB reference structure                                                                                                                                               |
| `ddG`                                                                                                                                                                                                                                                                                                                                 | Variation of the folding energy upon mutation calculated by FoldX (kcal/mol)                                                                                                                  |
| `SecStruc`                                                                                                                                                                                                                                                                                                                            | Secondary structure type in the reference structure                                                                                                                                           |
| `ddS`                                                                                                                                                                                                                                                                                                                                 | Variation of the vibrational entropy  upon mutation calculated from ENCom output (kcal/mol)                                                                                                   |
| `subScore`                                                                                                                                                                                                                                                                                                                            | BLOSUM62 substitution score                                                                                                                                                                   |
| `grantham`                                                                                                                                                                                                                                                                                                                            | Grantham score                                                                                                                                                                                |
| `sneath`                                                                                                                                                                                                                                                                                                                              | Sneath's index                                                                                                                                                                                |
| `classWT` and `classMut`                                                                                                                                                                                                                                                                                                              | Amino acid class of the reference and the variant respectively according to Albatineh and Razeghifard (2008). c1: M, V, I, L; c2 :G, A, C, P, S, T; c3 : F, W, Y; c4: R, H, K; c5: N, Q, D, E |
| `sasa_mut` and `sasa_wt`                                                                                                                                                                                                                                                                                                              | Solvent-accessible surface area of the reference and the variant structures respectively                                                                                                      |
| `hb_mut` and `hb_wt`                                                                                                                                                                                                                                                                                                                  | Number of hydrogen bonds established by the reference amino acid and the variant amino acid respectively                                                                                      |
| `sb_mut` and `sb_wt`                                                                                                                                                                                                                                                                                                                  | Number of salt bridges established by the reference amino acid and the variant amino acid respectively                                                                                        |
| `sasa_ratio`                                                                                                                                                                                                                                                                                                                          | Ratio of accessible surface area beween the reference and the variant amino acid.                                                                                                             |
| `hyrophob_WT` and `hyrophob_Mut`                                                                                                                                                                                                                                                                                                      | Hydrophobicity of the reference and the variant amino acids respectively according to JOND750101 descriptor from AAindex database                                                             |
| `volume_WT` and `volume_Mut`                                                                                                                                                                                                                                                                                                          |  Amino acid volume of the reference and variant residues respectively according to BIGC670101 descriptor from AAindex database                                                                |
| `pssm_wt` and `pssm_mut`                                                                                                                                                                                                                                                                                                              | PSSM scores of the reference and variant amino acids respectively.                                                                                                                            |
| `disulfide_breakage` `buried_Pro_introduced` `buried_glycine_replaced` `buried_hydrophilic_introduced` `buried_charge_introduced` `buried_charge_switch` `sec_struct_change` `buried_charge_replaced` `buried_exposed_switch` `exposed_hydrophilic_introduced` `Buried_salt_bridge_breakage` `Large_helical_penality_in_alpha_helix ` | Indicates if a significant structural change occured upon the mutation according to Missense3D criteria, 1: True, 0: False                                                                    |
| `swaat_prediction`                                                                                                                                                                                                                                                                                                                    | Predictive model outcome, 1: significant impact, 0: Neutral                                                                                                                                   |



## 3.2 HTML output

The HTML report is better in summarizing significant structural events with likely significant impact. These can be spotted as red flag tags. The report also links the genome coordinates with the amino acid substitution. In addition, the report offers the sequence annotation associated with each amino acid position as well as some statistics about the drug-binding  likelihood of the amino acid. These statistics were calculated from the FTMap output which reports the number of probe molecules that likely bind to the amino acid.  

# 4. Example: annotattion of six ADME genes

The example shows the annotation process of seven ADME genes  *CYP2A13* , *CYP2B6* , *CYP2C19*.  *CYP2D6*.  *CYP3A4*, and *CYP1A2*. We have optained the VCF files of these gene from [PharmVar](https://www.pharmvar.org/) and store them in `inputexample` directory. We generated a list of all the gene symbols corresponding to the VCF files in `./inputexample/gene_list.txt` shown below: 

```
CYP2A13
CYP2B6
CYP2C19
CYP2D6
CYP3A4
CYP1A2
```

The gene list will tell SWAAT to limit the annotation only for the gene listed in the file. 

In this case the user must specify the location of the VCF files, the location of the gene list file with an option to specify the path and name of the output folder. Here the final result is saved in `swaat_output`. The workflow can then be run from the command line. 

```bash
nextflow run main.nf --dbhome /home/hothman//SWAAT/database/ \
		--vcfhome /home/hothman/SWAAT/inputexample/VCFs \
		--outfolder ./swaat_out \
		--genelist ./inputexample/gene_list.txt
```

The workflow will process 113 variants in total. 

---

## Funding

This work was primarily funded through a grant by GlaxoSmithKline Research and Development Ltd to the Wits Health Consortium. The funders had no role in the design of the study; in the collection,  analyses, or interpretation of data; in the writing of the manuscript,  or in the decision to publish the results.

## Acknowledgement

We thank [Scott Hazelhurst](http://dept.ee.wits.ac.za/~scott/) and Nhlamulo Khoza from the University of the Witwatersrand who both spent significant effort to improve and test the tools developed in this project. 

## Authors 

Houcemeddine Othman, Sherlyn Jemimah and Jorge EB da Rocha.

Please contact Houcemeddine Othman ([houcemoo@gmail.com](mailto:houcemoo@gmail.com)/[houcemeddine.othman@wits.ac.za](mailto:houcemeddine.othman@wits.ac.za)) for any questions or comments.

