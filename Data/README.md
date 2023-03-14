# README for Cleaning the Raw Dataset
Date:2023/3/10  
Author: Guyuan Tang

### Description
The directory `/Data` contains the raw dataset `Eurasian_raw.xlsx` as well as the curated test dataset `test_Eurasian.xlsx` for testing our tool. The curation is finished by using the program `clean_data.py`.  

The `test_Eurasian.xlsx` contains only the needed columns for testing purpose, which are as the followings: Locality, Country, Lat.(latitude), Long.(longitude), Y_haplogroup, mt_group.


### Data
The raw data was named `Eurasian - Dataset_tims.xlsx`, which was downloaded from the BINP29 canvas page (https://canvas.education.lu.se/courses/22587/files/3514096?module_item_id=839243).  
It was renamed to `Eurasian_raw.xlsx` for convenience.


### Important Note!
To execute the program `clean_data.py`, the following packages/libraries should be installed to the environment in advance: `geopy` (v2.3.0) and `openpyxl` (v3.1.0).  
They could be installed via `conda`:  
(we suggest to create a new conda environment for better execution)
```shell
conda install geopy=2.3.0
conda install openpyxl=3.1.0
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
This program mainly reformatted some information in the raw data `Eurasian_raw.xlsx`:  
1) The missing values were represented as ".." in the raw data, and we corrected them into np.nan;  
2) Some individuals did not have information on longitude and latitude, we fulfilled the geography information using their location and country;
3) The formats for some columns (e.g. mt_haplogroup, country, Y_haplogroup) were not the same, we curated them;
4) The "Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]" was renamed into "Ages_2000" and we added 50 years to each to represent the date mean in BP in years before 2000CE. And we created a new column "Ages_interval" to group the individuals into 13 groups, with each group having a 1000 years period.
#### 3. Output for the program
The output for `clean_data.py` would be `test_Eurasian.xlsx` in `\Data` in the repository.  
It contains all the essential information for our tool `HaploMap.py` to execute. And the column names should be the same as `test_Eurasian.xlsx` for our tool to perform successfully.  
  
The essential information and the strict column names are as the following:  
"Ages_2000",
"Age_interval",
"Locality",  
"Country",  
"Lat.",  
"Long.",  
"Y_haplogroup",  
"mt_haplogroup"
