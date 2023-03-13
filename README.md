# README for HaploMap
Date: 2023/03/13  
Author: Guyuan Tang

## 1. Description
HaploMap is a tool for analyzing the haplogroups on Y-chromosome and mitochondrial DNA (mtDNA) on an input dataset, which could output meaningful information for the users. It has three modes as the followings:  
- mode 1: finding the closest haplogroup for the query and plotting the individuals (if any) belonging to both the query and closest haplogroups on the world map.  
- mode 2: searching for the haplogroup that was defined by the query mutation, providing mutation information, and plotting the individuals (if any) in the haplogroup on the map.  
- mode 3: calculating the main haplogroup frequencies in a selected country.  

### Folders and Files
The main script for our tool is `HaploMap.py` in the main page, along with the essential data file `SNP_index.xlsx` and a test dataset `test_Eurasian.xlsx`. A zip file `HaploMap.zip` is also provided for the users to download which includes all the essential files mentioned above.  

(1) The folder `\Data` contains: the raw dataset `Eurasian_raw.xlsx` downloaded from the BINP29 course webpage (https://canvas.education.lu.se/courses/22587/files/3514096?module_item_id=839243), a script `clean_data.py` to clean and to reformat the raw dataset, the example test dataset `test_Eurasian.xlsx` produced by the script and containing all the essential information as a suitable input for our tool.  
To note, `clean_data.py` was designed specifically for the raw dataset `Eurasian_raw.xlsx`, but it could be refered to when the users preparing their own dataset. See more information in the `README.md` within `\Data`.  

(2) The folder `\SNP_index` contains: the raw index file `SNP_Index_raw.xlsx` downloaded from the ISOGG website (https://docs.google.com/spreadsheets/d/1UY26FvLE3UmEmYFiXgOy0uezJi_wOut-V5TD0a_6-bE/edit#gid=1934392066), a script `clean_SNP.py` to remove missing values and to correct the format, the data `SNP_index_1.xlsx` produced by the script and being one of the essential database for our tool to refer to.  
To know more details, please read the `README.md` in the same folder.  

(3) The folders `Example_output_*` contain the example outputs from our tool in different mode, which used the same input file `test_Eurasian.xlsx`. Our tool could be executed in different working directories, just to make sure that the `HaploMap.py` and `SNP_index.xlsx` should be in the same directory.  
For example, the outputs in `Example_output_1` were produced when working in that directory.
```shell
python ../HaploMap.py --mode 1 --input ../Data/test_Eurasian.xlsx
```

## 2. Installation
To use HaploMap, just download the `HaploMap.zip` and uncompress it. Make sure to keep the `SNP_index.xlsx` in the same directory as `HaploMap.py`.  
### Dependencies
Before using the `HaploMap`, the users should install three packages in their working environment: `openpyxl`=3.1.0, `geopy`, `geopandas`.  
Here we provide an example for preparing the working environment by using `conda`.  
We suggest to create a new environment to avoid some unexpected conflicts.
```shell
# create a new environment in conda
conda create -n HaploMap
conda activate HaploMap
# install required packages, we also provide the versions that we have tested
conda install -c conda-forge openpyxl=3.1.0 # there is some bugs in the newer v3.1.1, and thus we recommend you to install v3.1.0!
conda install -c conda-forge geopy=0.12.2
conda install -c conda-forge geopands=0.12.2
```

## 3. Input Format
If the user do not want to use our test dataset, there are some reminders on the input format.  
The input for our tool must be an Excel file at this moment. The user could refer to the test dataset `test_Eurasian.xlsx` we provided. The column names are suggested to be the same as the test dataset.  
The input file should have at least the following information (column names are provided next to the information): country ("Country"), longitude ("Long."), latitude ("Lat."), haplogroup on Y chromosome ("Y_haplogroup"), haplogroup on mtDNA ("mt_haplogroup").  
Note: If your dataset only includes haplogroup on Y chromosome or on mtDNA, it is also fine. They do not need to be included in a same dataset. 


## 4. How to Use
The following examples on how to use our tool used the `test_Eurasian.xlsx` as the input.  
For help, type the command `HaploMap.py --help` in the terminal.

### (1) Mode 1
Go to the working directory where you want your outputs being stored. Here we use the `Example_output_1` as the example.  
```
cd Example_output_1
``` 

To use mode 1, simply type the command in the terminal:
```shell
python ../Haplomap --mode 1 --input ../Data/test_Eurasian.xlsx
```
And then follow the questions showing on the screen to input the query. There are two questions:  
- Please select the chromosome (Y/mt):  
- Enter the haplogroup:  

Note: in mode 1, our tool will search for the 3 top closest haplogroups of the query input by default. If there is not any individuals from the closest haplogroup in the dataset, it will search for the second closest (which is the closest one of the first closest haplogroup), and then the third closest. You could change the default value by altering the loop number `n` in the script.  

Note: The search will end in these two cases: a) the search reaches the origin haplogroup (Y-Adam chromosome for Y-DNA, mt_MRCA Eve DNA for mtDNA); b) there is not any individuals from the three (by default) top closest haplogroups in the dataset.
#### Output
The output graph would be named as: chromosome + the name of the closest haplogroup + ".pdf" in PDF format. For example, the `mt_W.pdf` would be an output with the query haplogroup W1 on mtDNA.  

The output would be a graph plotting the individuals from the query (if any) and the closest haplogroups. By default, the individuals are plotted in black and red, belonging to the query and the closest haplogroups, respectively. 

### (2) Mode 2
#### Output

### (3) Mode 3
#### Output

## 5. Limitations

## 6. FAQs