# README for HaploMap
Date: 2023/03/13  
Author: Guyuan Tang

## 1. Description
HaploMap is a tool for analyzing the haplogroups on Y-chromosome and mitochondrial DNA (mtDNA) on an input dataset, which could output meaningful information for the users. It has three modes as the followings:  
- mode 1: finding the closest haplogroup for the query and plotting the individuals (if any) belonging to both the query and closest haplogroups on the world map, seperated by the age intervals of the individuals.  
- mode 2: searching for the haplogroup that was defined by the query mutation, providing mutation information, and plotting the individuals (if any) in the haplogroup on the map.  
- mode 3: calculating the main haplogroup frequencies in a selected country.  

### Folders and Files
The main script for our tool is `HaploMap.py` in the main page, along with the essential data file `SNP_index.xlsx` and a test dataset `test_Eurasian.xlsx`. A zip file `HaploMap.zip` is also provided for the users to download which includes all the essential files mentioned above.  

(1) The folder `\Data` contains: the raw dataset `Eurasian_raw.xlsx` downloaded from the BINP29 course webpage (https://canvas.education.lu.se/courses/22587/files/3514096?module_item_id=839243), a script `clean_data.py` to clean and to reformat the raw dataset, the example test dataset `test_Eurasian.xlsx` produced by the script and containing all the essential information as a suitable input for our tool.  
To note, `clean_data.py` was designed specifically for the raw dataset `Eurasian_raw.xlsx`, but it could be refered to when the users preparing their own dataset. See more information in the `README.md` within `\Data`.  

(2) The folder `\SNP_index` contains: the raw index file `SNP_Index_raw.xlsx` downloaded from the ISOGG website (https://drive.google.com/drive/folders/1lnQ53cURAzTCOoyIUBIW3z1S-g41KrJc?usp=share_link), a script `clean_SNP.py` to remove missing values and to correct the format, the data `SNP_index_1.xlsx` produced by the script and being one of the essential database for our tool to refer to.  
To know more details, please read the `README.md` in the same folder.  

(3) The folders `Example_output_*` contain the example outputs from our tool in different mode, which used the same input file `test_Eurasian.xlsx`. Our tool could be executed in different working directories, just to make sure that the `HaploMap.py` and `SNP_index.xlsx` should be in the same directory.  
For example, the outputs in `Example_output_1` were produced when working in that directory.
```shell
python ../HaploMap.py --mode 1 --input ../Data/test_Eurasian.xlsx
```

(4) The folder `HaploMap_download` contains the `HaploMap.zip` (containing only `HaploMap.py` and `SNP_index.xlsx`) the `HaploMap-linux.zip` and `HaploMap-windows.zip` which were both packaged with all the dependencies the program needs. This folder will provide an external link (https://drive.google.com/file/d/13uPUkTEab7o9ewD0Kahg1Nx9uitlNNux/view?usp=sharing) to google drive for you to download the chosen file.

## 2. Installation
There are two ways of installing HaploMap: download the `HaploMap.zip` (includes the `HaploMap.py`, `SNP_index.xlsx`, `README.md`), or download the corresponding file from the `HaploMap_download` folder in the external drive based on your operating system.  
Here we provide guidance on these two ways:
### (1) Download the source code along with the index
We suggest you to use this method to install HaploMap if you want to use your own dataset and you could alter some settings in the source code.  
Or if you have already installed the dependencies (described below) in your working environment, you could use this method to avoid repeated installation of the packages.  
Users for MacOS could also use this method.

To use HaploMap, just download the `HaploMap.zip` and uncompress it. Make sure to keep the `SNP_index.xlsx` in the same directory as `HaploMap.py`.  
### Dependencies
Before using the `HaploMap`, the users should install two packages in their working environment: `openpyxl`=3.1.0 and `geopandas`.  
Here we provide an example for preparing the working environment by using `conda`.  
We suggest to create a new environment to avoid some unexpected conflicts.
```shell
# create a new environment in conda
conda create -n HaploMap
conda activate HaploMap
# install required packages, we also provide the versions that we have tested
conda install -c conda-forge openpyxl=3.1.0 # there is some bugs in the newer v3.1.1, and thus we recommend you to install v3.1.0!
conda install -c conda-forge geopands=0.12.2
```


### (2) Download the zip file based on the operating system (Linux or Windows)
We suggest you to use this method if you want to easily start using our tool.  
Go to external google drive: https://drive.google.com/file/d/13uPUkTEab7o9ewD0Kahg1Nx9uitlNNux/view?usp=sharing. Here you could find the folder containing the zip files.
Download the zip file in `HaploMap_download` based on your operating system.  
- for Windows user: download the `HaploMap-windows.zip`, uncompress it
- for Linux user: download the `HaploMap-linux.zip`, uncompress it  

Note: if you install our tool in this way, make sure to specify the directory where you download and uncompress the zip file. By default, the script will be stored in a folder called `HaploMap`.


## 3. Input Format
If the user do not want to use our test dataset, there are some reminders on the input format.  
The input for our tool must be an Excel file at this moment. The user could refer to the test dataset `test_Eurasian.xlsx` we provided. The column names are suggested to be the same as the test dataset.  
The input file should have at least the following information (column names are provided next to the information): country ("Country"), longitude ("Long."), latitude ("Lat."), haplogroup on Y chromosome ("Y_haplogroup"), haplogroup on mtDNA ("mt_haplogroup"), age-range for each individual ("Age_interval").  
Note: If your dataset only includes haplogroup on Y chromosome or on mtDNA, it is also fine. They do not need to be included in a same dataset.  


## 4. How to Use
The following examples on how to use our tool used the `test_Eurasian.xlsx` as the input.  
For help, type the command `HaploMap.py --help` in the terminal. 
When choosing the mode and specifying the input dataset, `--mode` command does not need to be before the `--input`. But make sure that `--mode` is followed by number 1, 2 or 3; while `--input` should be followed by the input dataset.  

Note: always be careful about the directory of the script, you need to specify it in the command line. You need to specify the directory where you uncompress the zip file, and add an additional `/HaploMap/HaploMap.py` in order to execute the program.

### (1) Mode 1
Go to the working directory where you want your outputs being stored. Here we use the `Example_output_1` as the example.  
```
cd Example_output_1
``` 

In this case, the script is located in `../`. The directory should be altered based on where you uncompress the zip file. For example, `../HaploMap/HaploMap.py`.  

To use mode 1, simply type the command in the terminal:
```shell
python ../HaploMap.py --mode 1 --input ../Data/test_Eurasian.xlsx
```
And then follow the questions showing on the screen to input the query. There are two questions:  
- Please select the chromosome (Y/mt):  
- Enter the haplogroup:  

Note: in mode 1, our tool will search for the 3 top closest haplogroups of the query input by default. If there is not any individuals from the closest haplogroup in the dataset, it will search for the second closest (which is the closest one of the first closest haplogroup), and then the third closest. You could change the default value by altering the loop number `n` in the script.  

Note: The search will end in these two cases: a) the search reaches the origin haplogroup (Y-Adam chromosome for Y-DNA, mt_MRCA Eve DNA for mtDNA); b) there is not any individuals from the three (by default) top closest haplogroups in the dataset.

#### Output
The output graph would be named as: chromosome + the name of the closest haplogroup + ".pdf" in PDF format. For example, the `mt_W.pdf` would be an output with the query haplogroup W1 on mtDNA.  

The output would be a graph plotting the individuals from the query (if any) and the closest haplogroups. By default, the individuals are plotted in triangles and circles, belonging to the query and the closest haplogroups, respectively. 

### (2) Mode 2  
Go to the working directory where you want your outputs being stored. Here we use the `Example_output_2` as the example.  
```
cd Example_output_2
``` 
In this case, the script is located in `../`. The directory should be altered based on where you uncompress the zip file. For example, `../HaploMap/HaploMap.py`.  

To use mode 2, simply type the command in the terminal:
```shell
python ../HaploMap.py --mode 2 --input ../Data/test_Eurasian.xlsx
```
And then follow the questions showing on the screen to input the query. There is one question:  
- Please enter the mutation name (Y-DNA):  

Note: always make sure that the file `SNP_index.xlsx` is in the same directory as the main program.

#### Output
The outputs will include a report on the mutation (including the haplogroup it defines, build 37 number, build 38 number, mutation information, and the number of individuals belonging to that haplogroup in the dataset).  

If there is not any individuals in the dataset found in the haplogroup defined by the input mutation, the output would only be a txt file reporting the mutation information.  
If there are individuals in the dataset found in the defined haplogroup, a graph similar to mode 1 will be printed to the working directory. It shows the location of the individuals on the map and seperates them by different colors according to their age intervals.

### (3) Mode 3
Go to the working directory where you want your outputs being stored. Here we use the `Example_output_3` as the example.  
```
cd Example_output_3
``` 
In this case, the script is located in `../`. The directory should be altered based on where you uncompress the zip file. For example, `../HaploMap/HaploMap.py`.  

To use mode 3, simply type the command in the terminal:
```shell
python ../HaploMap.py --mode 3 --input ../Data/test_Eurasian.xlsx
```
And then follow the questions showing on the screen to input the query. There are twp questions:  
- Please select the chromosome (Y/mt):
- Please select a country to discover:

Note: the input for country is case-sensitive, and thus be sure you have entered a correct country name. The abbreviations for countries are not accepted (for example, "UK" could not be recognized as "United Kingdom").

#### Output
The output for mode 3 would be a txt file reporting the haplogroup frequency within a country. Ages of the individuals are not considered here. And thus, this mode would be more meaningful in analyzing the haplogroup frequency on a populaiton within a short time period or the same time frame.

## 5. Limitations
At this moment, HaploMap has some limitations:
- Before using our tool, the users should be able to install two dependent packages which were used by our program. Or choose the second way to download the zip file with all the dependencies' files.
- MacOS users may need to choose the first method to download `HaploMap.zip` and prepare the environment (we suggest conda here). 
- The input file should be an excel file with the essential column names as our example dataset provided.  
- Some default settings in the script are designed specifically for the example dataset, although they could be altered in a simple way (for example, the colors for different age intervals).

## 6. Important Messages
(1) When using command line to execute the program, make sure the directory should be specified to where the `HaploMap.py` locates.  
(2) When not sure what to type in the command line to execute the program, take a careful look at the example above or type `python HaploMap.py --help` for help.  
(3) If you want to use your own dataset, please refer to the `clean_data.py` in the folder `/Data` to contain essential information with the right column names in your own dataset.  
(4) If you download the `HaploMap.zip` instead of the `HaploMap-linux.zip` or `HaploMap-windows.zip`, remember to keep the `HaploMap.py` together with the `SNP_index.xlsx` in the same directory.  
(5) If you have any questions, we welcome you to discuss with us to improve our tool. You could email to gu5747ta-s@student.lu.se.