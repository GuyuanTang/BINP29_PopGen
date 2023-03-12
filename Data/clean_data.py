# -*- coding: utf-8 -*-
#!/usr/bin/python3
"""
Title: clean_data.py
Date: 2023/03/08
Author: Guyuan Tang

Description: the code for preparing a test dataset suitable for the HaploMap.py to execute based on the original dataset Eurasian_tims.xlsx

List of packages/libraries:
    pandas, numpy, geopy, sys
    
    
List of functions:
    No self-created function.
    
    
Steps: 1) curate the missing values in Y and mtDNA columns; 2) in mt_DNA column, some individuals are classified as the "haplogroup + mutation points" which means their subclades are unclear, consider them as missing value because we want accurate data to find the closest groups; 3) correct the format of column "Country" (some may have one or more spaces at the end) and make curation to individuals lacking longitude and latitude; 4) export the dataframe to a new excel file.
    
    
    
Usage: python clean_data.py Eurasian_raw.xlsx test_Eurasian.xlsx
"""

import pandas as pd
import numpy as np
from geopy.geocoders import Nominatim
import sys

in_file = sys.argv[1]
out_file = sys.argv[2]


# read the raw data into a Dataframe, and extract the columns we need for this test dataset
df_raw = pd.read_excel(in_file,
                   usecols=["#", "Long.", "Lat.", "Locality", "Country", "Y haplogroup  in ISOGG v15.73 notation (automatically called)", "mtDNA haplogroup if ≥2 or published"],
                   header=0)

# rename the columns "Y haplogroup  in ISOGG v15.73 notation (automatically called)" and "mtDNA haplogroup if ≥2 or published" into easy-called column names: "Y_haplogroup" and "mt_haplogroup"
df = df_raw.rename(columns={'Y haplogroup  in ISOGG v15.73 notation (automatically called)':'Y_haplogroup', 'mtDNA haplogroup if ≥2 or published':'mt_haplogroup'})

# in the raw data, ".." was used to represent missing values in Y-haplogourp and mt_haplogroup
# replace the ".." into np.nan which could be recognized as missing values in python
df.loc[df['Y_haplogroup']=="..", ['Y_haplogroup']] = np.nan
df.loc[df['mt_haplogroup']=="..", ['mt_haplogroup']] = np.nan

# replace the the n/a into np.nan as missing values
df.loc[df['Y_haplogroup'].str.contains('n/a', regex=True, na=False), ['Y_haplogroup']] = np.nan
df.loc[df['mt_haplogroup'].str.contains('n/a', regex=True, na=False), ['mt_haplogroup']] = np.nan

# in mt_haplogroup column, individuals with unclear haplogroup (closest group + possible mutation sites) could be represented as missing values
df.loc[df['mt_haplogroup'].str.contains("[-\+\*/_\s]", regex=True, na=False), ['mt_haplogroup']] = np.nan
# further curate some rows with strange characters, such as U5a1a1¬†, M3a1b.., H1c5a'
df['mt_haplogroup'].replace(r"¬†","", regex=True, inplace=True)
df['mt_haplogroup'].replace(r"'$","", regex=True, inplace=True)
df['mt_haplogroup'].replace(r"\.\.","", regex=True, inplace=True)


# curate the format of Country, some may have a space at the end
df['Country'].replace(r"\s+$","", regex=True, inplace=True)

# curate the longitude and latitude for individuals lacking the information
geolocator = Nominatim(user_agent="MyApp")
# find the indice for individuals lacking longitude and latitude
noLLindex = df[df['Long.']=='..'].index.tolist()
for i in noLLindex:
    locality = df.loc[i,'Locality']
    country = df.loc[i,'Country']
    location = geolocator.geocode(locality)
    if location: # the locality could be recognized by the geopy
        df.loc[i,'Long.'] = location.longitude
        df.loc[i,'Lat.'] = location.latitude
    else: # the locality could not be recogized by the geopy and we use the country instead
        location = geolocator.geocode(country)
        df.loc[i,'Long.'] = location.longitude
        df.loc[i,'Lat.'] = location.latitude

# print the dataset into a new excel file
df.to_excel(out_file, index=False)
