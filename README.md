# SWAAT



## Preparing dependencies for the auxiliary workflow 

### Installing `tranvar`

```
# Use miniconda3 for a cleaner installation
pip install --user transvar

# set up databases
transvar config --download_anno --refversion hg19

# in case you don't have a reference
transvar config --download_ref --refversion hg19

# test
transvar panno -i 'PIK3CA:p.E545K' --ucsc --ccds
```

### Installing `BioPython`

```
pip install biopython
# Or if you prefer conda 
conda install -c anaconda biopython 
```

