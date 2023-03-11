# -*- coding: utf-8 -*-
#!/usr/bin/python3
"""
Title: clean_SNP.py
Date: 2023/03/11
Author: Guyuan Tang

Description: the code for preparing the SNP index on Y-DNA suitable for the HaploMap.py to execute

List of packages/libraries: pandas, sys
    
    
List of functions: No self-design function
    
Steps: 
    1) remove the empty rows (which are unclear in the index without build37 numbers);
    2) curate the mutation names (remove characters such as "^^");
    3) export to a new excel file.

Usage: python clean_SNP.py SNP_Index_raw.xlsx SNP_index_1.xlsx
"""

import pandas as pd
import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

df_raw = pd.read_excel(in_file, usecols=[0,1,4,5,6], skiprows=[0,2]) # the first row contains text information about the index, the third row contains unclear name of the mutation

# remove the rows without build37 number which are also unclear at the mutation names
df = df_raw.dropna(axis=0, subset=['Build 37 Number'])

# curate the mutation names (remove the "^^" characters)
df['Name'].replace(r"\^{1,2}$","", regex=True, inplace=True)

# export to a new excel file
df.to_excel(out_file, index=False)
