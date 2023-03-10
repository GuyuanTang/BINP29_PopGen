# README for Cleaning the Raw Dataset
Date:2023/3/10  
Author: Guyuan Tang

### Description
The directory `/Data` contains the raw dataset `Eurasian_raw.xlsx` as well as the curated test dataset `test_Eurasian.xlsx` for testing our tool. The curation is finished by using the program `clean_data.py`.  

The `test_Eurasian.xlsx` contains only the needed columns for testing purpose, which are as the followings: Locality, Country, Lat.(latitude), Long.(longitude), Y_haplogroup, mt_group.


### Data
The raw data was named `Eurasian - Dataset_tims.xlsx`, which was downloaded from the BINP29 canvas page (https://canvas.education.lu.se/courses/22587/files/3514096?module_item_id=839243).  
It was renamed to `Eurasian_raw.xlsx` for convenience.


### Important Note
To execute the program `clean_data.py`, the following packages/libraries should be installed to the environment in advance: `geopy` (v2.3.0) and `openpyxl` (v3.0.10).  
They could be installed via `conda`:  
(we suggest to create a new conda environment for better execution)
```shell
conda create -n popgen
conda activate popgen
conda install geopy=2.3.0
conda install openpyxl=3.0.10
```

The program was designed specifically for this `Eurasian_raw.xlsx` dataset, and it may not fit in with other datasets. But you could use it as a reference for preparing the dataset that are suitable as the input for our tool.


### Information on `clean_data.py`
#### 1. Usage
Open the working directory and type `python clean_data.py [input] [output]` in the terminal. The program only accepts excel file as the input.  

For example:
```
python clean_data.py Eurasian_raw.xlsx test_Eurasian.xlsx
```
#### 2. Script for the program
The script could be found in the `clean_data.py` in the same directory.
