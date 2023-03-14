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
    
    
Steps: 1) curate the missing values in Y and mtDNA columns; 2) in mt_DNA column, some individuals are classified as the "haplogroup + mutation points" which means their subclades are unclear, consider them as missing value because we want accurate data to find the closest groups; 3) calculating the individuals birth year ranges, which were seperated in to 13 groups from 0-11000 years (every 1000 years as a group) before 2000 CE; 4) correct the format of column "Country" (some may have one or more spaces at the end) and make curation to individuals lacking longitude and latitude; 5) export the dataframe to a new excel file.
    
    
    
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
                   usecols=["#", "Long.", "Lat.", "Locality", "Country", "Y haplogroup  in ISOGG v15.73 notation (automatically called)", "mtDNA haplogroup if ≥2 or published", "Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]"],
                   header=0)

# rename the columns "Y haplogroup  in ISOGG v15.73 notation (automatically called)" and "mtDNA haplogroup if ≥2 or published" into easy-called column names: "Y_haplogroup" and "mt_haplogroup"
df = df_raw.rename(columns={'Y haplogroup  in ISOGG v15.73 notation (automatically called)':'Y_haplogroup', 'mtDNA haplogroup if ≥2 or published':'mt_haplogroup', 'Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]':'Ages_2000'})

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

# Add 50 years to each "Ages_2000" to represent the sample's age from 2000, which was 1950 in the raw data
df['Ages_2000'] = df['Ages_2000'].map(lambda x: x+50)
# Create a new column to contain the approximate range for each individual's birth year
for year in range(0,10001,1000):
    year_range = year + 1000
    if (2000 - year > 0):  # calculating years in CE
        year_max = 2000 - year
        year_min = 2000 - year_range + 1
        interval = str(year_min) + "-" + str(year_max) + " CE"
        df.loc[((df['Ages_2000']>=year) & (df['Ages_2000']<year_range)), ['Age_interval']] = interval
    
    else: # calculating years in BCE
        year_max = abs(2000 - year) + 1
        year_min = abs(2000 - year_range) # there is no 0 cE or 0 BCE, it goes from 1 BCE to 1 CE in the calender
        interval = str(year_min) + "-" + str(year_max) +" BCE"
        df.loc[((df['Ages_2000']>=year) & (df['Ages_2000']<year_range)), ['Age_interval']] = interval
        

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
