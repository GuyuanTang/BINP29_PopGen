# README for Cleaning the Raw Dataset
Date:2023/3/11  
Author: Guyuan Tang

### Description
The directory `/SNP_index` contains the raw dataset `SNP_Index_raw.xlsx` as well as the curated dataset `SNP_index_1.xlsx`. The curation is finished by using the program `clean_SNP.py`.  

The `SNP_index_1.xlsx` contains only the needed columns for testing purpose, which are as the followings: Name, Subgroup Name, Build37 Number, Build 38 Number, Mutation Info.


### Data
The raw data was named `SNP Index .xlsx`, which was downloaded from the ISOGG webpage (https://docs.google.com/spreadsheets/d/1UY26FvLE3UmEmYFiXgOy0uezJi_wOut-V5TD0a_6-bE/edit#gid=1934392066).  
Reference: International Society of Genetic Genealogy. Y-DNA Haplogroup Tree 2019, Version: 15.73, Date: 11 July 2020, http://www.isogg.org/tree/ [Date of access: 9 Mar 2023]  

It was renamed to `SNP_Index_raw.xlsx` for convenience.


### Important Note
To execute the program `clean_SNP.py`, the following package should be installed to the environment in advance: `openpyxl` (v3.1.0).  
They could be installed via `conda`:  
```shell
conda activate popgen
conda install openpyxl=3.1.0
```

The program was designed specifically for this `SNP_Index_raw.xlsx` dataset, and it may not fit in with other datasets. But you could use it as a reference for preparing the dataset that are suitable as the input for our tool.


### Information on `clean_SNP.py`
#### 1. Usage
Open the working directory and type `python clean_SNP.py [input] [output]` in the terminal. The program only accepts excel file as the input.  

For example:
```
python clean_SNP.py SNP_Index_raw.xlsx SNP_index_1.xlsx
```
A warning may occur but it will not affect our desired output.
#### 2. Script for the program
The script could be found in the `clean_SNP.py` in the same directory.  
The program removed the mutation that were unclear (most without the GRCh37 - Build 37 Number). And it also curated the mutation names that contained the "^^". 