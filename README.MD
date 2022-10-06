# PROT TOOLS
## A collection of tools to map crosslinks with jwalk
### 2019-09 Kai Kammer

## Installation
### Versions
The scripts were tested to be working with the following package versions:

* python: 3.10.5
* pandas: 1.4.3
* numpy: 1.23.1
* biopyton: 1.79
* pkai: 1.1
* jwalk: 1.3 according to the python file name or 1.0.0 according to the setup.py

### Conda Environment

It is recommended to install all required packages into their own separate environment. Conda/Mamba are recommended but a pip venv should work fine as well.

    conda create -n prot_tools
    conda activate prot_tools
    conda install ipython numpy pandas biopython setuptools_scm pytorch -c conda-forge
    pip install pkai

* the models need to be downloaded by hand
* find your site-packages path with

    pip show pkai | grep Location  
    wget https://github.com/bayer-science-for-a-better-life/pKAI/raw/main/pKAI/models/pKAI_model.pt -P your-path/site-packages/pKAI/models/  
    wget https://github.com/bayer-science-for-a-better-life/pKAI/raw/main/pKAI/models/pKAI+_model.pt -P your-path/site-packages/pKAI/models/


## Usage

### Fasta PDB Mapping
First map your protein's fasta file to its pdb file with

    prot_match_pdb_chains.py

For example:

    prot_match_pdb_chains.py my_prot.pdb my_prot.fasta

Multiple fasta and pdb files are possible which are validated via biopyton. Note that offsets between pdb and fasta are computed automatically on a per-chain basis. The above command will create a file called prot_uni_to_pdb_chain.csv.

### pKa Calculation
pKa values are calculated with pKAI (https://github.com/bayer-science-for-a-better-life/pKAI). 
The script

    prot_get_pka.py

is a small wrapper for pKAI to create output compatible with the following steps. It also support the offsets file created in the previous step:

    prot_get_pka.py my_prot.pdb -off prot_uni_to_pdb_chain.csv

This will create a file call prot_pka.csv.

### SASA
Solvent-accessible surface area is calculated via biopython's PDB.SASA.ShrakeRupley module.
The script named

    prot_get_sasa.py

formats the output into a form compatible with the following steps. It also support the offsets file created in the previous step:

    prot_get_sasa.py my_prot.pdb -off prot_uni_to_pdb_chain.csv

This will create a file call prot_sasa.csv.

### SASD
Solvent-accessible surface distance is calculated via jwalk. Either use the online server https://jwalk.ismb.lon.ac.uk/jwalk/ or create a separate conda environment (as jwalk requires python 2).

    conda create -n jwalk
    conda activate jwalk
    conda install numpy biopython python=2.7.18
    cd /tmp
    git clone https://github.com/Topf-Lab/Jwalk
    cd Jwalk/Jwalk
    pip install .

Jwalk can then be used like this:

    jwalk -i my_prot.pdb -max_dist 35

where max_dist limits the euclidean maximum distance between two amino acids (which is always shorter than the SASD).  
Jwalk will create a Jwalk_results folder. The file we need will be called my_prot_crosslink_list.txt.  
The output written to that file then formatted and accorded for fasta/pdb offsets with

    prot_jwalk_format_output.py

For example

    prot_jwalk_format_output.py Jwalk_results/my_prot_crosslink_list.txt -off prot_uni_to_pdb_chain.csv

will write to a file called prot_distances.csv

## Summary
After the steps above you should end up with 3 files:

* prot_pka.csv
* prot_sasa.csv
* prot_distances.csv

These files are ready to use as input for creating a kinetic model (via xy.py)
