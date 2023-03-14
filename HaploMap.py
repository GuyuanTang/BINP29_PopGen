# -*- coding: utf-8 -*-
#!/usr/bin/python3
"""
Title: HaploMap.py
Date: 2023/03/10
Author: Guyuan Tang

Description: the code for analysing the population genetics on human Y-chromosome and mtDNA. It mainly has three functions/modes: 
    1) select a haplogroup to plot the individuals within it and its closest haplogroup (if there is no individuals within the cloest haplogroup, it will search for the second closest and then the third closest) on the map; 
    2) select a mutation name (SNP) on Y-chromosome to get its haplogroup and relevant information (bulid 37,38, and mutation) and plot the individuals (if any) on the map; 
    3) select a country and chromosome to report its haplogroup frequency (represented as the main haplogroups instead of their subclades).

Reference for constructing the phylogenic trees of the haplogroups:
    1) Y-DNA: International Society of Genetic Genealogy. Y-DNA Haplogroup Tree 2019, Version: 15.73, Date: 11 July 2020, http://www.isogg.org/tree/ [Date of access: 9 Mar 2023].
    2) mt-DNA: van Oven M, Kayser M. 2009. Updated comprehensive phylogenetic tree of global human mitochondrial DNA variation. Hum Mutat 30(2):E386-E394. http://www.phylotree.org. DOI:10.1002/humu.20921. [Date of access: 9 Mar 2023].
    
List of packages/libraries: 
    sys, pandas, numpy, geopandas, matplotlib, re
    
List of functions:
    get_close(dict_data, haplogroup)
    map_plot(out_file, haplo_df, query_df, match_df, query, haplo_close)
    haplo_freq(df, chr_name, country_name)

Steps are described in each mode.    

Example Usage: python HaploMap.py --mode [1/2/3] --input test_Eurasian.xlsx
-Mode 1: input a haplogroup and search for the closest ones in the dataset and plot them on the map (it will search for 3 times at maximum)
And follow the instruction on the screen to enter the chromosome (Y/mt) and the haplogroup name.

-Mode 2: input a mutation name on Y-chromosome to get its haplogroup and relevant information, and plot the individuals (if any) on the map
And follow the instruction on the screen to enter the mutation name.

-Mode 3: input a country and chromosome to report its haplogroup frequency
And follow the instruction on the screen to enter the country name.

For help: python HaploMap.py --help

Note: the column names are suggested to be set as the test dataset (test_Eurasian.xlsx). Please refer to it for detailed format.
    
"""

import pandas as pd
import numpy as np
import geopandas as gpd
import re
import matplotlib.pyplot as plt
import sys

if len(sys.argv) == 2 and sys.argv[1] == '--help':
    print("\nExample Usage: python HaploMap.py --mode [1/2/3] --input test_Eurasian.xlsx\n\n-Mode 1: input a haplogroup and search for the closest ones in the dataset and plot them on the map (it will search for 3 times at maximum)\nAnd follow the instruction on the screen to enter the chromosome (Y/mt) and the haplogroup name.\n\n-Mode 2: input a mutation name on Y-chromosome to get its haplogroup and relevant information, and plot the individuals (if any) on the map\nAnd follow the instruction on the screen to enter the mutation name.\n\n-Mode 3: input a country and chromosome to report its haplogroup frequency\nAnd follow the instruction on the screen to enter the country name.")

elif len(sys.argv) == 5 and ('--mode' in sys.argv) and ('--input' in sys.argv):
    mode_loc = sys.argv.index('--mode')
    input_loc = sys.argv.index('--input')
    mode = sys.argv[mode_loc+1]
    in_file = sys.argv[input_loc+1]
    

    """
1. The first and the main function: input a haplogroup, search for its closest ones and plot them on the map.

Steps:
    1) develop a function get_close(dict_data, haplogroup) to search for the dictionary to get the closest group
    2) develop a function map_plot(out_file, query_df, match_df, query, haplo_close) for plotting on the map
    
    -for Y-DNA:
        1) Create a dictionary to contain the main tree trunk as well as a list on the haplogroups that do not follow the regular naming rule (for a haplogroup to develop its subclade, a number or a lower-case letter was added to the end, whether it was a number or a letter depended on the end character);
        2) If the input is an element in the main tree, it will call the get_close() to find its closest haplogroup; if it is not in the main tree, it will find the closest haplogroup following the naming rule (but we set both the confirm and approximate(ending with ~) groups as our searching target, and only one of them would be found in the dataset);
        3) If individuals were found, we extracted their indice to create a dataframe for plotting; if there is not any individual, the program will search for the second closest haplogroup and etc. (we set a maximum value at 3 by default, the variable "num").
    
    -for mt-DNA:
        1) Create dictionaries to contain both the main tree and all the subgroups. Also, a list was also created as to link the subgroups to the main tree trunk;
        2) The get_close() function would be called to find the closest haplogroups;
        3) If individuals were found, we extracted their indice to create a dataframe for plotting; if there is not any individual, the program will search for the second closest haplogroup and etc. (we set a maximum value at 3 by default, the variable "num").


Example usage: 
    -for Y-DNA: python HaploMap.py --mode 1 --input test_Eurasian.xlsx 
        Please select the chromosome(Y/mt): Y 
        Enter the haplogroup: R1a1a1b1a1
    -for mt_DNA: python HaploMap.py --mode 1 --input test_Eurasian.xlsx
        Please select the chromosome(Y/mt): mt 
        Enter the haplogroup: U5a1a1a

        """

    if mode == '1':
        chr_name = input("Please select the chromosome (Y/mt):")
        haplo_name = input("Enter the haplogroup:")
        num = 3
        
        
        # define a color list for plotting
        intervals = ['1001-2000 CE', '1-1000 CE', '1000-1 BCE', '2000-1001 BCE', '3000-2001 BCE', '4000-3001 BCE', '5000-4001 BCE', '6000-5001 BCE', '7000-6001 BCE', '8000-7001 BCE', '9000-8001 BCE', '10000-9001 BCE', '11000-10001 BCE']
        color_list = ['lightcoral','brown','red','darkorange','gold','yellowgreen','limegreen','blue','violet','fuchsia','darkorchid','yellow','cyan']
        colors = {}
        for i in range(len(intervals)):
            key = intervals[i]
            colors[key] = color_list[i]
        
        # define a function to get the closest haplogroup
        def get_close(dict_data, haplogroup):
            for k,v in dict_data.items():
                if haplogroup in v:
                    return k
        
        # define a function for plotting on map
        def map_plot(out_file, query_df, match_df, query, haplo_close, colors):
            # start plotting
            world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
            fig,ax=plt.subplots(dpi=600)
            ax.set_aspect('equal')
            ax.xaxis.label.set_visible(False)
            ax.yaxis.label.set_visible(False)
            world.plot(ax=ax, color = 'lightgrey', edgecolor = 'darkgrey', linewidth = 0.5)
            if not query_df.empty: # only plot the query individuals if there is any
                grouped_query = query_df.groupby('Age_interval')
                for key, group in grouped_query:
                    label_name = query + " (" + key + ")"
                    group.plot(ax=ax, kind='scatter', x='Long.', y='Lat.', s=3.5, label=label_name, color=colors[key], marker = '^')
            # plot the individuals from the closest haplogroup
            grouped = match_df.groupby('Age_interval')
            for key, group in grouped:
                label_name = haplo_close + " (" + key + ")"
                group.plot(ax=ax, kind='scatter', x='Long.', y='Lat.', s=3.5, label=label_name, color=colors[key])
            plt.legend(bbox_to_anchor=(1.0,1.0), loc=2, fontsize=8)
            plt.savefig(out_file, format="pdf", bbox_inches='tight')
        
        
        # For Y-DNA
        if chr_name.lower() == 'y':
    
            # Create a dictionary to contain main tree trunk on Y-DNA
            # the haplogroups on Y-DNA have a naming rule to follow, and thus we only store the information on the main tree trunk
            # for the haplogroup K, we contain the patterns from other haplogroups (e.g. L,T,NO,etc.) instead of calling them K1,K2a,etc..
            Y_main = ['A0000','A000-T','A000','A00-T','A00','A0-T','A0','A1','A1b','A1b1','BT','B','CT','DE','CF','D','E','C','F','GHIJK','G','HIJK','H','IJK','IJ','K','I','J','LT','K2','L','T','NO','K2b','NO1','N','O','K2b1','P','S','M','P1','Q','R']
            Y_tree = {'Y':['A0000','A000-T'],
          'A000-T':['A000','A00-T'], 'A00-T':['A00','A0-T'], 'A0-T':['A0','A1'], 'A1':['A1b'], 'A1b':['A1b1','BT'],
          'BT':['B','CT'],
          'CT':['DE','CF'],
          'DE':['D','E'],
          'CF':['C','F'],
          'F':['GHIJK'],
          'GHIJK':['G','HIJK'],
          'HIJK':['H','IJK'],
          'IJK':['IJ','K'],
          'IJ':['I','J'],
          'K':['LT','K2'], 'LT':['L','T'], 'K2':['NO','K2b'], 'NO':['NO1'], 'NO1':['N','O'], 'K2b':['K2b1','P'], 'K2b1':['S','M'],
          'P':['P1'], 'P1':['Q','R']}
            
            try:
                # read the input Eurasian dataset
                haplo_df = pd.read_excel(in_file, header=0)
                
                # repeat the input from the user
                print("\nYou have selected {} on chromosome Y".format(haplo_name))
                # search for the individuals belong to the query haplogroup
                query_index = haplo_df[haplo_df['Y_haplogroup']==haplo_name].index.tolist()
                # create a dataframe to contain the selected individuals for later plotting
                query_df = haplo_df.iloc[query_index]
                # check if there is any individual belong to the query haplogroup
                if not query_index:
                    print("\nThere isn't any individual belong to {} in the dataset.".format(haplo_name))
                query = haplo_name
                    
                # by default, the program will only search for three times (if no match in the first closest, it will seach for the second and the third)
                n = 0
                while n<num:
                    # first check if the input haplogroup is in the main trunk (Y-main)
                    if haplo_name in Y_main:
                        haplo_close = get_close(dict_data=Y_tree, haplogroup=haplo_name)
                        print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        if haplo_close == "Y":
                            n = num
                            print("\nThe search reached the end: Adam Y chromosome! Please try another haplogroup.")
                        else:
                            # create a list to contain the indice for the matched individuals
                            individual = haplo_df[haplo_df["Y_haplogroup"]==haplo_close].index.tolist()
                        
                            # if there are individuals belong to the closest haplogroup
                            if individual:
                                n = num
                                match_df = haplo_df.iloc[individual]
                                individual_count = len(match_df)
                                print("Found {} individual(s) in the haplogroup {}".format(individual_count,haplo_close))
                                
                                # create the output file name
                                out_file = "Y_" + haplo_close + ".pdf"
                                
                                # plot the results on the map
                                map_plot(out_file, query_df, match_df, query, haplo_close, colors)
                                print("\nThank you for using HaploMap! The graph has printed to your working directory.")
                                
                            else: # there isn't any individual in the closest group
                                n += 1
                                if n<num:
                                    print("There isn't any individual in the closest haplogroup, searching for next closest!")    
                                    haplo_name = haplo_close
                                else:
                                    print('\nThank you for using HaploMap! We only search for the top {} closest haplogroups. No matched individuals!'.format(num))
                        
                    # if the haplo_name is not in the main trunk
                    else:
                        # according to the rule for naming the haplogroup, we remove the final character/number for each search
                        if haplo_name.endswith("~"):
                            # in this case, we need to remove two characters
                            # some haplogroups only have approximate locations in the tree, and thus we also include them in the search
                            haplo_close_0 = haplo_name[:-2] # with confirmed location
                            haplo_close_1 = haplo_close_0 + "~" # with approximate location
                            if (haplo_df['Y_haplogroup'].eq(haplo_close_0)).any() or (haplo_df['Y_haplogroup'].eq(haplo_close_1)).any():
                                # only one haplo_close_* would match
                                n = num
                                if (haplo_df['Y_haplogroup'].eq(haplo_close_0)).any():
                                    print("The closest group {} is found!".format(haplo_close_0))
                                    individual = haplo_df[haplo_df["Y_haplogroup"]==haplo_close_0].index.tolist()
                                    match_df = haplo_df.iloc[individual]
                                    individual_count = len(match_df)
                                    print("Found {} individual(s) in the haplogroup {}".format(individual_count,haplo_close_0))
                                    # create the output file name
                                    out_file = "Y_" + haplo_close_0 + ".pdf"
                                    # plot the results on the map
                                    map_plot(out_file, query_df, match_df, query, haplo_close_0, colors)
                                    print("\nThank you for using HaploMap! The graph has been printed to your working directory.")
                                
                                elif (haplo_df['Y_haplogroup'].eq(haplo_close_1)).any():
                                    print("The closest group {} is found!".format(haplo_close_1))
                                    individual = haplo_df[haplo_df["Y_haplogroup"]==haplo_close_1].index.tolist()
                                    match_df = haplo_df.iloc[individual]
                                    individual_count = len(match_df)
                                    print("Found {} individual(s) in the haplogroup {}".format(individual_count,haplo_close_1))
                                    # create the output file name
                                    out_file = "Y_" + haplo_close_1 + ".pdf"
                                    # plot the results on the map
                                    map_plot(out_file, query_df, match_df, query, haplo_close_1, colors)
                                    print("\nThank you for using HaploMap! The graph has been printed to your working directory.")
                            # do not have match to the closest haplogroup
                            else:
                                n += 1
                                if n<num:
                                    print("There isn't any individual in the closest haplogroup, searching for next closest!")    
                                    haplo_name = haplo_close_0
                                else:
                                    print('\nThank you for using HaploMap! We only search for the top {} closest haplogroups. No matched individuals!'.format(num))
                                
                            
                        else:
                            haplo_close_0 = haplo_name[:-1] # with confirmed location
                            haplo_close_1 = haplo_close_0 + "~" # with approximate location
                            if (haplo_df['Y_haplogroup'].eq(haplo_close_0)).any() or (haplo_df['Y_haplogroup'].eq(haplo_close_1)).any():
                                # only one haplo_close_* would match
                                n = num
                                if (haplo_df['Y_haplogroup'].eq(haplo_close_0)).any():
                                    print("The closest haplogroup {} is found!".format(haplo_close_0))
                                    individual = haplo_df[haplo_df["Y_haplogroup"]==haplo_close_0].index.tolist()
                                    match_df = haplo_df.iloc[individual]
                                    individual_count = len(match_df)
                                    print("Found {} individual(s) in the haplogroup {}".format(individual_count,haplo_close_0))
                                    # create the output file name
                                    out_file = "Y_" + haplo_close_0 + ".pdf"
                                    # plot the results on the map
                                    map_plot(out_file, query_df, match_df, query, haplo_close_0, colors)
                                    print("\nThank you for using HaploMap! The graph has been printed to your working directory.")
                                
                                elif (haplo_df['Y_haplogroup'].eq(haplo_close_1)).any():
                                    print("The closest haplogroup {} is found!".format(haplo_close_1))
                                    individual = haplo_df[haplo_df["Y_haplogroup"]==haplo_close_1].index.tolist()
                                    match_df = haplo_df.iloc[individual]
                                    individual_count = len(match_df)
                                    print("Found {} individual(s) in the haplogroup {}".format(individual_count,haplo_close_1))
                                    # create the output file name
                                    out_file = "Y_" + haplo_close_1 + ".pdf"
                                    # plot the results on the map
                                    map_plot(out_file, query_df, match_df, query, haplo_close_1, colors)
                                    print("\nThank you for using HaploMap! The graph has been printed to your working directory.")
                            # do not have match to the closest haplogroup
                            else:
                                n += 1
                                if n<num:
                                    print("There isn't any individual in the closest haplogroup, searching for next closest!")    
                                    haplo_name = haplo_close_0
                                else:
                                    print('\nThank you for using HaploMap! We only search for the top {} closest haplogroups. No matched individuals!'.format(num))
            
            except FileNotFoundError as not_found:
                print("The file {} was not found!".format(not_found.filename))
            
            except IndexError:
                print("Please check the haplogroup input again!")
        
        
        
        
        # For mtDNA
        elif chr_name.lower() == 'mt':
            # Create a dictionary to contain full tree trunk on mtDNA
            # the naming rule for defining haplogroups on mtDNA is different in some subclades, and thus we store the information on the whole tree trunk
            mt_main = ["L1'2'3'4'5'6","L0","L2'3'4'5'6","L1","L2'3'4'6","L5","L3'4'6","L2","L3'4","L6","L3","L4","N","M","Q","CZ","C","Z","E","G","D","R","A","O","S","X","I","W","Y","P","U","J","T","HV","H","V","B4'5","B6","F","JT","K"]

            mt_tree_main = {'mt-MRCA':["L1'2'3'4'5'6","L0"],
           "L1'2'3'4'5'6":["L2'3'4'5'6","L1"],
           "L2'3'4'5'6":["L2'3'4'6","L5"],
           "L2'3'4'6":["L3'4'6","L2"],
           "L3'4'6":["L3'4","L6"],
           "L3'4":["L3","L4"],
           "L3":["N","M"], 
           "M":["M8","M9","M12'G","M29'Q","M80'D"],
           "M29'Q":["Q"], 
           "M8":["CZ"], "CZ":["C","Z"],
           "M9":["E"],
           "M12'G":["G"],
           "M80'D":["D"],
           "N":["R","N1'5","N2","N9","A","O","S","X"],
           "N1'5":["N1"],"N1":["N1a"], "N1a":["N1a1'2"], "N1a1'2":["N1a1"], "N1a1":["N1a1b"],"N1a1b":["I"],
           "N2":["W"],
           "N9":["Y"],
           "R":["R0","R2'JT","R9","R11'B6","B4'5","P","U"],
           "R2'JT":["R2","JT"], "JT":["J","T"],
           "R11'B6":["B6"],
           "R0":["HV"], "HV":["HV0","H"], "HV0":["HV0a"], "HV0a":["V"],
           "R9":["F"],
           "U":["U2'3'4'7'8'9"],"U2'3'4'7'8'9":["U8"], "U8":["U8b'c"], "U8b'c":["U8b"], "U8b":["K"]
           }

            mt_tree_L = {'L0':["L0a'b'f'g'k", "L0d"], "L0a'b'f'g'k":["L0a'b'f'g", "L0k"], "L0a'b'f'g":["L0a'b'g","L0f"], "L0a'b'g":["L0a'g","L0b"], "L0a'g":["L0a","L0g"], "L0a":["L0a1'4","L0a2","L0a3"], "L0a1'4":["L0a1","L0a4"], "L0a1":["L0a1a","L0a1b","L0a1c","L0a1d","L0a1e"], "L0a1a":["L0a1a1","L0a1a2","L0a1a3"], "L0a1b":["L0a1b1","L0a1b2"], "L0a1b1":["L0a1b1a"], "L0a1b1a":["L0a1b1a1"], "L0a1b1a1":["L0a1b1a1a"], "L0a1b2":["L0a1b2a"], "L0a1c":["L0a1c1"], "L0a2":["L0a2a","L0a2b","L0a2c","L0a2d"], "L0a2a":["L0a2a1","L0a2a2"], "L0a2a1":["L0a2a1a","L0a2a1b"], "L0a2a1a":["L0a2a1a1","L0a2a1a2"], "L0a2a2":["L0a2a2a"], "L0a2a2a":["L0a2a2a1","L0a2a2a2"], "L0a2b":["L0a2b1"], "L0f":["L0f1","L0f2"], "L0f2":["L0f2a","L0f2b"], "L0f2a":["L0f2a1"], "L0k":["l0k1","L0k2"], "L0k1":["L0k1a","L0k1b"], "L0k1a":["L0k1a1","L0k1a2","L0k1a3"], "L0k1a1":["L0k1a1a","L0k1a1b","L0k1a1c","L0k1a1d"], "L0k1a2":["L0k1a2a"], "L0k2":["L0k2a","L0k2b"], "L0k2a":["L0k2a1"], "L0k2a1":["L0k2a1a"], "L0d":["L0d1'2","L0d3"], "L0d1'2":["L0d1","L0d2"], "L0d1":["L0d1a'c'd","L0d1b"], "L0d1a'c'd":["L0d1a'd","L0d1c"], "L0d1a'd":["L0d1a","L0d1d"], "L0d1a":["L0d1a1"], "L0d1a1":["L0d1a1a","L0d1a1b","L0d1a1c","L0d1a1d"], "L0d1a1a":["L0d1a1a1","L0d1a1a2","L0d1a1a3"], "L0d1a1b":["L0d1a1b1"], "L0d1a1b1":["L0d1a1b1a","L0d1a1b1b"], "L0d1c":["L0d1c1","L0d1c2","L0d1c3"], "L0d1c1":["L0d1c1a"], "L0d1c1a":["L0d1c1a1","L0d1c1a2"], "L0d1c1a1":["L0d1c1a1a","L0d1c1a1b"], "L0d1c1a1a":["L0d1c1a1a1","L0d1c1a1a2"], "L0d1c2":["L0d1c2a"], "L0d1c2a":["L0d1c2a1"], "L0d1b":["L0d1b1","L0d1b2"], "L0d1b1":["L0d1b1a","L0d1b1b","L0d1b1c"], "L0d1b1a":["L0d1b1a1"], "L0d1b1b":["L0d1b1b1"], "L0d1b2":["L0d1b2a","L0d1b2b"], "L0d1b2a":["L0d1b2a1","L0d1b2a2"], "L0d1b2b":["L0d1b2b1","L0d1b2b2"], "L0d1b2b1":["L0d1b2b1a","L0d1b2b1b"], "L0d1b2b1b":["L0d1b2b1b1"], "L0d1b2b2":["L0d1b2b2a","L0d1b2b2b","L0d1b2b2c"], "L0d1b2b2b":["L0d1b2b2b1"], "L0d1b2b2c":["L0d11b2b2c1","L0d1b2b2c2"], "L0d2":["L0d2a'b'd","L0d2c"], "L0d2a'b'd":["L0d2a","L0d2b","L0d2d"], "L0d2a":["L0d2a1","L0d2a2"], "L0d2a1":["L0d2a1a","L0d2a1b","L0d2a1c"], "L0d2a1a":["L0d2a1a1","L0d2a1a2","L0d2a1a3"], "L0d2a1a1":["L0d2a1a1a"], "L0d2b":["L0d2b1","L0d2b2"], "L0d2b1":["L0d2b1a","L0d2b1b"], "L0d2b1a":["L0d2b1a1"], "L0d2b1a1":["L0d2b1a1a"], "L0d2c":["L0d2c1","L0d2c2"], "L0d2c1":["L0d2c1a","L0d2c1b"], "L0d2c1a":["L0d2c1a1"], "L0d2c2":["L0d2c2a","L0d2c2b"], "L0d2c2a":["L0d2c2a1"], "L0d2c2a1":["L0d2c2a1a"], "L0d3":["L0d3a","L0d3b"], "L0d3b":["L0d3b1","L0d3b2"],
           "L1":["L1b","l1c"], "L1b":["L1b1","L1b2'3"], "L1b1":["L1b1a"], "L1b1a":["L1b1a1'4","L1b1a2","L1b1a3","L1b1a9","L1b1a15","L1b1a17","L1b1a18","L1b1a5","L1b1a6","L1b1a7","L1b1a8","L1b1a10","L1b1a12","L1b1a13","L1b1a14","L1b1a16"], "L1b1a1'4":["L1b1a1","L1b1a4"], "L1b1a4":["L1b1a4a"], "L1b1a2":["L1b1a2a"], "L1b1a3":["L1b1a3a","L1b1a3b"], "L1b1a3a":["L1b1a3a1"], "L1b1a15":["L1b1a15a"], "L1b1a7":["L1b1a7a"], "L1b1a10":["L1b1a10a","L1b1a10b"], "L1b1a12":["L1b1a12a","L1b1a12b"], "L1b2'3":["L1b2","L1b3"], "L1b2":["L1b2a"], "L1c":["L1c1'2'4'5'6","L1c3"], "L1c1'2'4'5'6":["L1c1'2'4'6","L1c5"], "L1c1'2'4'6":["L1c1","L1c2'4","L1c6"], "L1c1":["L1c1a'b'd","L1c1c"], "L1c1a'b'd":["L1c1a","L1c1b'd"], "L1c1a":["L1c1a1","L1c1a2"], "L1c1a1":["L1c1a1a","L1c1a1b"], "L1c1a1a":["L1c1a1a1","L1c1a1a2"], "L1c1a1a1":["L1c1a1a1a","L1c1a1a1b"], "L1c1a1a1b":["L1c1a1a1b1"], "L1c1a2":["L1c1a2a","L1c1a2b","L1c1a2c"], "L1c1a2a":["L1c1a2a1","L1c1a2a2"], "L1c1b'd":["L1c1b","L1c1d"], "L1c1b":["L1c1b1"], "L1c1d":["L1c1d1"], "L1c2'4":["L1c2", "L1c4"], "L1c2":["L1c2a","L1c2b"], "L1c2a":["L1c2a1","L1c2a2","L1c2a3"], "L1c2a1":["L1c2a1a","L1c2a1b"], "L1c2a3":["L1c2a3a"], "L1c2b":["L1c2b1","L1c2b2"], "L1c2b1":["L1c2b1a'b","L1c2b1c"], "L1c2b1a'b":["L1c2b1a","L1c2b1b"], "L1c2b1a":["L1c2b1a1"], "L1c2b1b":["L1c2b1b1"], "L1c4":["L1c4a","L1c4b"], "L1c3":["L1c3a","L1c3b'c"], "L1c3a":["L1c3a1"], "L1c3a1":["L1c3a1a","L1c3a1b"], "L1c3b'd":["L1c3b","L1c3c"], "L1c3b":["L1c3b1","L1c3b2"], "L1c3b1":["L1c3b1a","L1c3b1b"],
           "L5":["L5a","L5b"], "L5a":["L5a1","L5a2"], "L5a1":["L5a1a","L5a1b","L5a1c"], "L5b":["L5b1","l5b2"], "L5b1":["L5b1a","L5b1b"],
           "L2":["L2a'b'c'd","L2e"], "L2a'b'c'd":["L2a","L2b'c'd"], "L2a":["L2a1'2'3'4","L2a5"], "L2a1'2'3'4":["L2a1","L2a2'3'4"], "L2a1":["L2a1a","L2a1b","L2a1f","L2a1g","L2a1c","L2a1d","L2a1h","L2a1e","L2a1i","L2a1q","L2a1j","L2a1k","L2a1l","L2a1m","L2a1n","L2a1o","L2a1p"], "L2a1c":["L2a1c1","L2a1c6","L2a1c2","L2a1c3","L2a1c4","L2a1c5"], "L2a1c1":["L2a1c1a"], "L2a1c1a":["L2a1c1a1","L2a1c1a2"], "L2a1c2":["L2a1c2a"], "L2a1c3":["L2a1c3a","L2a1c3b"], "L2a1c3a":["L2a1c3a1"], "L2a1c3b":["L2a1c3b1","L2a1c3b2"], "L2a1c4":["L2a1c4a"], "L2a1c4a":["L2a1c4a1"], "L2a1d":["L2a1d1","L2a1d2"], "L2a1e":["L2a1e1"], "L2a1i":["L2a1i1"], "L2a1l":["L2a1l1","L2a1l2","L2a1l3"], "L2a1l1":["L2a1l1a","L2a1l1b"], "L2a1l1a":["L2a1l1a1"], "L2a1l2":["L2a1l2a"], "L2a1l2a":["L2a1l2a1"], "L2a1m":["L2a1m1"], "L2a1m1":["L2a1m1a"], "L2a2'3'4":["L2a2'3", "L2a4"], "L2a2'3":["L2a2","L2a3"], "L2a2":["L2a2a","L2a2b"], "L2a2a":["L2a2a1"], "L2a2b":["L2a2b1","L2a2b2"], "L2a2b1":["L2a2b1a"], "L2a4":["L2a4a","L2a4b"], "L2b'c'd":["L2b'c","L2d"], "L2b'c":["L2b","L2c"], "L2b":["L2b1","L2b2","L2b3"], "L2b1":["L2b1a","L2b1b"], "L2b1a":["L2b1a2","L2b1a3","L2b1a4"], "L2b2":["L2b2a"], "L2b3":["L2b3a","L2b3b","L2b3c"], "L2c":["L2c1","L2c2","L2c3","L2c4","L2c5"], "L2c1":["L2c1a"], "L2c2":["L2c2a","L2c2b"], "L2c2a":["L2c2a1"], "L2c2b":["L2c2b1","L2c2b2"], "L2c2b1":["L2c2b1a","L2c2b1b"], "L2c3":["L2c3a"], "L2d":["L2d1"], "L2d1":["L2d1a"], "L2e":["L2e1"], "L2e1":["L2e1a"],
           "L6":["L6a","L6b"],
           "L4":["L4a","L4b"], "L4a":["L4a1","L4a2"], "L4a1":["L4a1a"], "L4b":["L4b1","L4b2"], "L4b1":["L4b1a"], "L4b2":["L4b2a","L4b2b"], "L4b2a":["L4b2a1","L4b2a2"], "L4b2a2":["L4b2a2a","L4b2a2b","L4b2a2c"], "L4b2b":["L4b2b1"],
           "L3":["N","M","L3a","L3b'f","L3c'd","L3e'i'k'x","L3h"], "L3a":["L3a1","L3a2"], "L3a1":["L3a1a","L3a1b"], "L3a2":["L3a2a"], "L3b'f":["L3b","L3f"], "L3b":["L3b1","L3b2","L3b3"], "L3b1":["L3b1a","L3b1b"], "L3b1a":["L3b1a1","L3b1a2", "L3b1a3","L3b1a4","L3b1a5","L3b1a6","L3b1a7","L3b1a8","L3b1a9","L3b1a10","L3b1a11"], "L3b1a1":["L3b1a1a"], "L3b1a5":["L3b1a5a"], "L3b1a7":["L3b1a7a"], "L3b1a9":["L3b1a9a"], "L3b1b":["L3b1b1"], "L3b2":["L3b2a","L3b2b"], "L3f":["L3f1", "L3f2", "L3f3"], "L3f1":["L3f1a", "L3f1b"], "L3f1a":["L3f1a1"], "L3f1b":["L3f1b6","L3f1b1","L3f1b2","L3f1b3","L3f1b4","L3f1b5"], "L3f1b1":["L3f1b1a"], "L3f1b1a":["L3f1b1a1"], "L3f1b2":["L3f1b2a"], "L3f1b4":["L3f1b4a","L3f1b4b","L3f1b4c"], "L3f1b4a":["L3f1b4a1"], "L3f2":["L3f2a","L3f2b"], "L3f2a":["L3f2a1"], "L3f2a1":["L3f2a1a"], "L3f3":["L3f3a","L3f3b"], "L3c'd":["L3c","L3d"], "L3d":["L3d1'2'3'4'5'6"], "L3d1'2'3'4'5'6":["L3d1", "L3d2","L3d3","L3d4","L3d5","L3d6"], "L3d1":["L3d1a","L3d1b","L3d1c","L3d1d"], "L3d1a":["L3d1a1'2"], "L3d1a1'2":["L3d1a1","L3d1a2"], "L3d1a1":["L3d1a1a","L3d1a1b"], "L3d1a1a":["L3d1a1a1"], "L3d1b":["L3d1b1","L3d1b2","L3d1b3"], "L3d1b1":["L3d1b1a","L3d1b1b"], "L3d1b3":["L3d1b3a"], "L3d1c":["L3d1c1"], "L3d2":["L3d2a","L3d2b"], "L3d3":["L3d3a","L3d3b"], "L3d3a":["L3d3a1"], "L3d3a1":["L3d3a1a","L3d3a1b"], "L3d4":["L3d4a"], "L3d5":["L3d5a"], "L3e'i'k'x":["L3e","L3i","L3k","L3x"], "L3e":["L3e1","L3e2","L3e3'4'5"], "L3e1":["L3e1a","L3e1b","L3e1c","L3e1d","L3e1e","L3e1f","L3e1g"], "L3e1a":["L3e1a1","L3e1a2","L3e1a3"], "L3e1a1":["L3e1a1a"], "L3e1a3":["L3e1a3a","L3e1a3b"], "L3e1b":["L3e1b1","L3e1b2"], "L3e1d":["L3e1d1"], "L3e1d1":["L3e1d1a"], "L3e1e":["L3e1e1","L3e1e2"], "L3e1f":["L3e1f1","L3e1f2"], "L3e1f1":["L3e1f1a"], "L3e2":["L3e2a","L3e2b"], "L3e2a":["L3e2a1","L3e2a2","L3e2a3"], "L3e2a1":["L3e2a1a","L3e2a2b"], "L3e2a1b":["L3e2a1b1","L3e2a1b2","L3e2a1b3"], "L3e2b":["L3e2b1","L3e2b2","L3e2b3","L3e2b4","L3e2b5","L3e2b6","L3e2b7","L3e2b8"], "L3e2b1":["L3e2b1a"], "L3e2b1a":["L3e2b1a1","L3e2b1a2"], "L3e3'4'5":["L3e3'4","L3e5"], "L3e3'4":["L3e3","L3e4"], "L3e3":["L3e3a","L3e3b"], "L3e3b":["L3e3b1","L3e3b2","L3e3b3"], "L3e4":["L3e4a"], "L3e4a":["L3e4a1"], "L3e5":["L3e5a"], "L3e5a":["L3e5a1"], "L3e5a1":["L3e5a1a"], "L3i":["L3i1","L3i2"], "L3i1":["L3i1a","L3i1b"], "L3k":["L3k1"], "L3x":["L3x1","L3x2"], "L3x1":["L3x1a","L3x1b"], "L3x1a":["L3x1a1","L3x1a2"], "L3x2":["L3x2a","L3x2b"], "L3x2a":["L3x2a1"], "L3x2a1":["L3x2a1a"], "L3h":["L3h1","L3h2"], "L3h1":["L3h1a","L3h1b"], "L3h1a":["L3h1a1","L3h1a2"], "L3h1a2":["L3h1a2a","L3h1a3b"], "L3h1a2a":["L3h1a2a1"], "L3h1b":["L3h1b1","L3h1b2"], "L3h1b1":["L3h1b1a"]}

            mt_tree_M = {"M":["M1'20'51","M2","M3","M4''67","M5","M6","M7","M8","M9","M10","M11","M12'G", "M13'46'61","M14","M15","M17","M19'53","M21","M22","M23'75", "24'41","M25","M26","M27","M28","M29'Q","M31","M32'56","M33","M34'57","M35","M36","M39'70","M40","M42'74","M44","M47","M48","M49","M50","M52","M55'77", "M58","M59","M60","M62'68","M69","M71","M72","M73'79","M76","M81","M91","M80'D"],
"M1'20'51":["M1","M20","M51"], "M1":["M1a","M1b"], "M1a":["M1a1","M1a2","M1a3","M1a4","M1a5","M1a6","M1a7","M1a8"], "M1a1":["M1a1a","M1a1b","M1a1c","M1a1d","M1a1e","M1a1f","M1a1g","M1a1h","M1a1i"], "M1a1a":["M1a1a1"], "M1a1b":["M1a1b1","M1a1b2"], "M1a1b1":["M1a1b1a","M1a1b1b","M1a1b1c"], "M1a1b1b":["M1a1b1b1"], "M1a1e":["M1a1e1","M1a1e2"], "M1a2":["M1a2a","M1a2b"], "M1a3":["M1a3a","M1a3b"], "M1a3b":["M1a3b1","M1a3b2"], "M1a4":["M1a4a"], "M1a8":["M1a8a"], "M1b":["M1b1","M1b2"], "M1b1":["M1b1a","M1b1b"], "M1b2":["M1b2a","M1b2b","M1b2c"], "M51":["M51a","M51b"], "M51a":["M51a1","M51a2"], "M51a1":["M51a1a","M51a1b"], "M51b":["M51b1"], "M51b1":["M51b1a","M51b1b"], "M2":["M2a'b","M2c"], "M2a'b":["M2a","M2b"], "M2a":["M2a1","M2a2","M2a3"], "M2a1":["M2a1a","M2a1b","M2a1c"], "M2a1a":["M2a1a1","M2a1a2","M2a1a3"], "M2a1a1":["M2a1a1a","M2a1a1b"], "M2a1a1a":["M2a1a1a1"], "M2a1a1b":["M2a1a1b1"], "M2a1a2":["M2a1a2a"], "M2a1a2a":["M2a1a2a1"], "M2a1a2a1":["M2a1a2a1a"], "M2a1a3":["M2a1a3a","M2a1a3b"], "M2a1a3a":["M2a1a3a1"], "M2a2":["M2a2a"], "M2a3":["M2a3a"], "M2b":["M2b1","M2b2","M2b3","M2b4"], "M2b1":["M2b1a","M2b1b"], "M2b3":["M2b3a"], "M3":["M3a","M3b","M3c","M3d"], "M3a":["M3a1","M3a2"], "M3a1":["M3a1a","M3a1b"], "M3a2":["M3a2a"], "M3c":["M3c1","M3c2"], "M3c1":["M3c1a","M3c1b"], "M3c1b":["M3c1b1"], "M3c1b1":["M3c1b1a","M3c1b1b"], "M3d":["M3d1"],"M3d1":["M3d1a"], "M3d1a":["M3d1a1"], "M4''67":["M4","M65","M67", "M18'38", "M30", "M37","M43","M45","M54","M63","M64","M66"], "M4":["M4a","M4b"], "M65":["M65a","M65b"], "M65a":["M65a1","M65a2"], "M18'38":["M18","M38"], "M18":["M18a","M18b","M18c"], "M38":["M38a","M38b","M38c","M38d","M38e"], "M30":["M30a","M30b","M30c","M30d","M30e","M30f","M30g"], "M30a":["M30a1","M30a2"], "M30c":["M30c1"], "M30c1":["M30c1a"], "M30c1a":["M30c1a1"], "M30d":["M30d1","M30d2"], "M37":["M37a", "M37d","M37e"], "M37a":["M37a1"], "M37e":["M37e2"], "M43":["M43a","M43b"], "M43a":["M43a1"], "M45":["M45a"], "M66":["M66a", "M66b"], "M5":["M5a'd","M5b'c"], "M5a'd":["M5a","M5d"], "M5a":["M5a1","M5a2","M5a3","M5a4","M5a5"], "M5a1":["M5a1a","M5a1b"], "M5a2":["M5a2a"], "M5a2a":["M5a2a1","M5a2a2","M5a2a3", "M5a2a4"], "M5a2a1":["M5a2a1a"], "M5a2a1a":["M5a2a1a1","M5a2a1a2"], "M5a3":["M5a3a","M5a3b"], "M5b'c":["M5b","M5c"], "M5b":["M5b1","M5b2"], "M5b2":["M5b2a","M5b2b"], "M5b2b":["M5b2b1"], "M5b2b1":["M5b2b1a"], "M5c":["M5c1","M5c2"], "M6":["M6a","M6b"], "M6a":["M6a1","M6a2"], "M6a1":["M6a1a","M6a1b"], "M10":["M10a"], "M10a":["M10a1","M10a2"], "M10a1":["M10a1a","M10a1b"], "M10a1a":["M10a1a1"], "M10a1a1":["M10a1a1a","M10a1a1b"], "M10a1a1b":["M10a1a1b1","M10a1a1b2"], "M11":["M11a'b","M11c","M11d"], "M11a'b":["M11a","M11b"], "M11a":["M11a1","M11a2"], "M11b":["M11b1","M11b2"], "M11b1":["M11b1a"], "M11b1a":["M11b1a1"], "M12'G":["M12","G"], "M12":["M12a", "M12b"], "M12a":["M12a1","M12a2"], "M12a1":["M12a1a","M12a1b"], "M12a1a":["M12a1a1","M12a1a2"], "M12b":["M12b1","M12b2"], "M12b1":["M12b1a","M12b1b"], "M12b1a":["M12b1a1","M12b1a2"], "M12b1a2":["M12b1a2a","M12b1a2b"], "M12b2":["M12b2a"], "M13'46'61":["M13","M46","M61"], "M13":["M13a'b","M13c"], "M13a'b":["M13a","M13b"], "M13a":["M13a1","M13a2"], "M13a1":["M13a1a","M13a1b"], "M13a1b":["M13a1b1"], "M13b":["M13b1","M13b2"], "M46":["M46a"], "M61":["M61a"], "M17":["M17a","M17c"], "M17c":["M17c1"], "M17c1":["M17c1a"], "M17c1a":["M17c1a1"], "M17c1a1":["M17c1a1a"], "M19'53":["M19","M53"], "M53":["M53b"], "M21":["M21a","M21b"], "M21b":["M21b1","M21b2"], "M21b1":["M21b1a"], "M22":["M22a","M22b"], "M23'75":["M23","M75"], "M24'41":["M24","M41"], "M24":["M24a","M24b"], "M41":["M41a","M41b","M41c"], "M41a":["M41a1"], "M27":["M27a","M27b","M27c"], "M27a":["M27a1","M27a2","M27a3"], "M27a1":["M27a1a","M27a1b"], "M27a1a":["M27a1a1","M27a1a2"], "M27a2":["M27a2a","M27a2b"], "M27b":["M27b1","M27b2"], "M27b2":["M27b2a","M27b2b","M27b2c"], "M27b2a":["M27b2a1"], "M27b2b":["M27b2b1"], "M28":["M28a","M28b"], "M28a":["M28a1", "M28a2","M28a3","M28a4","M28a5","M28a6","M28a7"], "M28a2":["M28a2a"], "M28a5":["M28a5a","M28a5b"], "M28a6":["M28a6a"], "M28a7":["M28a7a","M28a7b"], "M28b":["M28b1"], "M29'Q":["M29","Q"], "M29":["M29a","M29b"], "M29b":["M29b1"], "M31":["M31a","M31b'c"], "M31a":["M31a1","M31a2"], "M31a1":["M31a1a","M31a1b"], "M31b'c":["M31b","M31c"], "M31b":["M31b1","M31b2"], "M32'56":["M32","M56"], "M32":["M32a","M32c"], "M33":["M33a","M33b","M33c","M33d"], "M33a":["M33a1","M33a2'3"], "M33a1":["M33a1a","M33a1b"], "M33a2'3":["M33a2","M33a3"], "M33a2":["M33a2a"], "M33a3":["M33a3a"], "M33b":["M33b1","M33b2"], "M34'57":["M34","M57"], "M34":["M34a","M34b"], "M34a":["M34a1","M34a2"], "M34a1":["M34a1a"], "M57":["M57a","M57b"], "M57b":["M57b1"], "M35":["M35a","M35b","M35c"], "M35a":["M35a1","M35a2"], "M35a1":["M35a1a"], "M35b":["M35b1","M35b2","M35b3","M35b4"], "M36":["M36a","M36b","M36c","M36d"], "M36d":["M36d1"], "M39'70":["M39","M70"], "M39":["M39a","M39b","M29c"], "M39a":["M39a1","M39a2"], "M39b":["M39b1","M39b2"], "M40":["M40a"], "M40a":["M40a1"], "M40a1":["M40a1a","M30a1b"], "M42'74":["M42","M74"], "M42":["M42a","M42b"], "M42b":["M42b1","M42b2"], "M42b1":["M42b1a"], "M74":["M74a","M74b"], "M74b":["M74b1","M74b2"], "M44":["M44a"], "M44a":["M44a1"], "M49":["M49a","M49c","M49d","M49e"], "M49a":["M49a1","M49a2"], "M49c":["M49c1"], "M49e":["M49e1"], "M50":["M50a"], "M50a":["M50a1","M50a2"], "M52":["M52a","M52b"], "M52a":["M52a1"], "M52a1":["M52a1a","M52a1b"], "M52a1b":["M52a1b1"], "M52b":["M52b1"], "M52b1":["M52b1a"], "M55'77":["M55","M77"], "M60":["M60a","M60b"], "M60a":["M60a1","M60a2"], "M62'68":["M62","M68"], "M62":["M62a","M62b"], "M62b":["M62b1","M62b2"], "M62b1":["M62b1a"], "M62b1a":["M62b1a1"], "M68":["M68a"], "M68a":["M68a1","M68a2"], "M68a1":["M68a1a"], "M68a2":["M68a2a"], "M69":["M69a"], "M71":["M71a","M71b","M71c"], "M71a":["M71a1","M71a2"], "M71a1":["M71a1a"], "M72":["M72a"], "M73'79":["M73","M79"], "M73":["M73a","M73b"], "M73a":["M73a1"], "M76":["M76a"], "M91":["M91a","M91b"], "M80'D":["M80","D"],
"M7":["M7a","M7b'c"], "M7a":["M7a1","M7a2"], "M7a1":["M7a1a","M7a1b"], "M7a1a":["M7a1a1","M7a1a2","M7a1a3","M7a1a4","M7a1a5","M7a1a6","M7a1a7","M7a1a8","M7a1a9"], "M7a1a1":["M7a1a1a"], "M7a1a4":["M7a1a4a"], "M7a1a5":["M7a1a5a"], "M7a1a6":["M7a1a6a"], "M7a1b":["M7a1b1","M7a1b2"], "M7a2":["M7a2a"], "M7a2a":["M7a2a1","M7a2a2","M7a2a3"], "M7a2a3":["M7a2a3a"], "M7b'c":["M7b","M7c"], "M7b":["M7b1","M7b2"], "M7b1":["M7b1a","M7b1b"], "M7b1a":["M7b1a1","M7b1a2"], "M7b1a1":["M7b1a1a","M7b1a1b","M7b1a1c","M7b1a1d","M7b1a1e","M7b1a1f","M7b1a1g","M7b1a1h","M7b1a1i"], "M7b1a1a":["M7b1a1a1","M7b1a1a2","M7b1a1a3"], "M7b1a1a1":["M7b1a1a1a","M7b1a1a1b","M7b1a1a1c","M7b1a1a1d"], "M7b1a1a1b":["M7b1a1a1b1"], "M7b1a1c":["M7b1a1c1"], "M7b1a1d":["M7b1a1d1"], "M7b1a1e":["M7b1a1e1","M7b1a1e2"], "M7b1a1i":["M7b1a1i1"], "M7b1a2":["M7b1a2a"], "M7b1a2a":["M7b1a2a1"], "M7b1a2a1":["M7b1a2a1a","M7b1a2a1b"], "M7b1a2a1b":["M7b1a2a1b1"], "M7b2":["M7b2a"], "M7c":["M7c1","M7c2","M7c3"], "M7c1":["M7c1a","M7c1b","M7c1c"], "M7c1a":["M7c1a1","M7c1a2","M7c1a3","M7c1a4","M7c1a5"], "M7c1a1":["M7c1a1a","M7c1a1b"], "M7c1a1a":["M7c1a1a1"], "M7c1a1b":["M7c1a1b1"], "M7c1a2":["M7c1a2a"], "M7c1a2a":["M7c1a2a1"], "M7c1a3":["M7c1a3a"], "M7c1a4":["M7c1a4a","M7c1a4b"], "M7c1b":["M7c1b1","M7c1b2"], "M7c1b2":["M7c1b2a","M7c1b2b"], "M7c1c":["M7c1c1","M7c1c2","M7c1c3"], "M7c1c1":["M7c1c1a"], "M7c1c1a":["M7c1c1a1"], "M7c1c2":["M7c1c2a"], "M7c1c3":["M7c1c3a","M7c1c3b","M7c1c3c","M7c1c3d","M7c1c3e","M7c1c3f","M7c1c3g","M7c1c3h","M7c1c3i"], "M7c1c3a":["M7c1c3a1"], "M7c2":["M7c2a","M7c2b"],
"M8":["M8a","CZ"], "M8a":["M8a1","M8a2'3"], "M8a1":["M8a1a"], "M8a2'3":["M8a2","M8a3"], "M8a2":["M8a2a","M8a2b","M8a2c","M8a2d","M8a2e"], "M8a2a":["M8a2a1"], "M8a3":["M8a3a"], "M8a3a":["M8a3a1"],
"M9":["M9a'b","E"], "M9a'b":["M9a","M9b"], "M9a":["M9a1","M9a4","M9a5"], "M9a1":["M9a1a","M9a1b"], "M9a1a":["M9a1a1","M9a1a2","M9a13"], "M9a1a1":["M9a1a1a","M9a1a1b","M9a1a1c","M9a1a1d"], "M9a1a1c":["M9a1a1c1"], "M9a1a1c1":["M9a1a1c1a","M9a1a1c1b","M9a1a1c1c"], "M9a1a1c1b":["M9a1a1c1b1","M9a1a1c1b2"], "M9a1a1c1b1":["M9a1a1c1b1a"], "M9a1a1c1b1a":["M9a1a1c1b1a1","M9a1a1c1b1a2"], "M9a1b":["M9a1b1","M9a1b2"], "M9a1b1":["M9a1b1a","M9a1b1b","M9a1b1c"], "M9a1b1a":["M9a1b1a1"], "M9a4":["M9a4a","M9a4b"], "M9a4a":["M9a4a1","M9a4a2"]}

            mt_tree_Q = {"Q":["Q1'2","Q3"], "Q1'2":["Q1","Q2"], "Q1":["Q1a","Q1b","Q1c","Q1d","Q1e","Q1f"], "Q1a":["Q1a1"], "Q1a1":["Q1a1a"], "Q1c":["Q1c1","Q1c2"], "Q1c1":["Q1c1a"], "Q1c2":["Q1c2a"], "Q1e":["Q1e1"], "Q1e1":["Q1e1a","Q1e1b","Q1e1c"], "Q1e1a":["Q1e1a1"], "Q1e1b":["Q1e1b1"], "Q1f":["Q1f1","Q1f2"], "Q2":["Q2a","Q2b"], "Q2a":["Q2a1","Q2a2","Q2a3","Q2a4"], "Q2a2":["Q2a2a","Q2a2b"], "Q2a3":["Q2a3a","Q2a3b"], "Q3":["Q3a","Q3b"], "Q3a":["Q3a1"]}

            mt_tree_C = {"C":["C1","C4","C5","C7"], "C1":["C1a","C1b","C1c","C1d","C1e","C1f"], "C1b":["C1b1","C1b2","C1b3","C1b4","C1b5","C1b6","C1b7","C1b10","C1b8","C1b9","C1b11","C1b12","C1b13","C1b14"], "C1b5":["C1b5a","C1b5b"], "C1b7":["C1b7a"], "C1b8":["C1b8a"], "C1b13":["C1b13a","C1b13b","C1b13c","C1b13d","C1b13e"], "C1b13a":["C1b13a1"], "C1b13c":["C1b13c1"], "C1c":["C1c1","C1c2","C1c3","C1c4","C1c5","C1c6","C1c7","C1c8"], "C1c1":["C1c1a","C1c1b"], "C1d":["C1d1","C1d2","C1d3"], "C1d1":["C1d1a","C1d1b","C1d1c","C1d1d"], "C1d1a":["C1d1a1"], "C1d1b":["C1d1b1"], "C1d1c":["C1d1c1"], "C1d2":["C1d2a"], "C4":["C4a'b'c", "C4d","C4e"], "C4a'b'c":["C4a","C4b","C4c"], "C4a":["C4a1","C4a2"], "C4a1":["C4a1a","C4a1b"], "C4a1a":["C4a1a1","C4a1a2","C4a1a3","C4a1a4","C4a1a5","C4a1a6"], "C4a1a1":["C4a1a1a"], "C4a1a2":["C4a1a2a"], "C4a1a3":["C4a1a3a","C4a1a3b","C4a1a3c","C4a1a3d"], "C4a1a3a":["C4a1a3a1"], "C4a1a4":["C4a1a4a"], "C4a2":["C4a2a","C4a2b","C4a2c"], "C4a2a":["C4a2a1"], "C4a2a1":["C4a2a1a","C4a2a1b"], "C4a2b":["C4a2b1","C4a2b2"], "C4a2b2":["C4a2b2a"], "C4a2c":["C4a2c1","C4a2c2"], "C4a2c2":["C4a2c2a"], "C4b":["C4b1","C4b2","C4b3","C4b5","C4b6","C4b7","C4b8"], "C4b1":["C4b1a","C4b1b"], "C4b2":["C4b2a"], "C4b3":["C4b3a", "C4b3b"], "C4b3a":["C4b3a1"], "C4b8":["C4b8a"], "C4c":["C4c1","C4c2"], "C4c1":["C4c1a","C4c1b"], "C5":["C5a","C5b","C5c","C5d"], "C5a":["C5a1","C5a2"], "C5a2":["C5a2a","C5a2b"], "C5a2b":["C5a2b1"], "C5b":["C5b1"], "C5b1":["C5b1a","C5b1b"], "C5b1a":["C5b1a1"], "C5b1b":["C5b1b1"], "C5c":["C5c1"], "C5c1":["C5c1a"], "C5d":["C5d1","C5d2"], "C7":["C7a","C7b"], "C7a":["C7a1","C7a2"], "C7a1":["C7a1a","C7a1c","C7a1d"], "C7a1a":["C7a1a1","C7a1a2"], "C7a2":["C7a2a"]}

            mt_tree_Z = {"Z":["Z1","Z2","Z3","Z4","Z5","Z7"], "Z1":["Z1a"], "Z1a":["Z1a1","Z1a2","Z1a3"], "Z1a1":["Z1a1a","Z1a1b"], "Z1a2":["Z1a2a"], "Z3":["Z3a","Z3b","Z3c","Z3d"], "Z3a":["Z3a1","Z3a2"], "Z3a1":["Z3a1a"], "Z4":["Z4a"], "Z4a":["Z4a1"], "Z4a1":["Z4a1a"], "Z4a1a":["Z4a1a1"]}

            mt_tree_E = {"E":["E1","E2"], "E1":["E1a"], "E1a":["E1a1","E1a2"], "E1a1":["E1a1a","E1a1b","E1a1c"], "E1a1a":["E1a1a1"], "E1a1a1":["E1a1a1a","E1a1a1b","E1a1a1c"], "E1a1a1b":["E1a1a1b1","E1a1a1b2"], "E1a1b":["E1a1b1","E1a1b2","E1a1b3","E1a1b4"], "E1a2":["E1a2a"], "E1a2a":["E1a2a1","E1a2a2","E1a2a3","E1a2a4"], "E2":["E2a","E2b"], "E2a":["E2a1","E2a2"], "E2a1":["E2a1a"], "E2b":["E2b1","E2b2"]}

            mt_tree_G = {"G":["G1","G2","G3","G4"], "G1":["G1a","G1b","G1c"], "G1a":["G1a1","G1a2'3"], "G1a1":["G1a1a","G1a1b"], "G1a1a":["G1a1a1","G1a1a2","G1a1a3","G1a1a4"], "G1a2'3":["G1a2","G1a3"], "G1b":["G1b1","G1b2","G1b3","G1b4"], "G1c":["G1c1","G1c2"], "G2":["G2a'c","G2b"], "G2a'c":["G2a","G2c"], "G2a":["G2a1","G2a2","G2a3","G2a4","G2a5"], "G2a1":["G2a1b","G2a1c","G2a1d","G2a1e","G2a1f","G2a1g","G2a1h"], "G2a1c":["G2a1c1","G2a1c2"], "G2a1d":["G2a1d1","G2a1d2"], "G2a1d1":["G2a1d1a"], "G2a1d2":["G2a1d2a"], "G2a1f":["G2a1f1"], "G2a2":["G2a2a"], "G2a3":["G2a3a"], "G2b":["G2b1","G2b2"], "G2b1":["G2b1a","G2b1b"], "G2b1a":["G2b1a1","G2b1a2"], "G2b2":["G2b2a","G2b2b","G2b2c"], "G3":["G3a","G3b"], "G3a":["G3a1'2","G3a3"], "G3a1'2":["G3a1","G3a2"], "G3a1":["G3a1a"], "G3a2":["G3a2a"], "G3b":["G3b1","G3b2"]}

            mt_tree_D = {"D":["D4","D5","D6"], "D4":["D1","D4a","D4b","D4c","D4d","D4e","D4f","D4g","D4h","D4i","D4j","D4k","D4o","D4p","D4l","D4m","D4n","D4q","D4s","D4t"], "D1":["D1a","D1b","D1c","D1d","D1e","D1f","D1g","D1h","D1i","D1j","D1k","D1m","D1n"], "D1a":["D1a1","D1a2"], "D1d":["D1d1","D1d2"], "D1f":["D1f1","D1f2","D1f3"], "D1g":["D1g1","D1g2","D1g5","D1g3","D1g4","D1g6"], "D1g1":["D1g1a","D1g1b"], "D1g2":["D1g2a"], "D1h":["D1h1","D1h2"], "D1i":["D1i1","D1i2"], "D1j":["D1j1"], "D1j1":["D1j1a"], "D1j1a":["D1j1a1","D1j1a2"], "D4a":["D4a1","D4a2","D4a3","D4a4","D4a5","D4a6","D4a7","D4a8"], "D4a1":["D4a1a","D4a1b","D4a1c","D4a1d","D4a1e","D4a1f","D4a1g","D4a1h"], "D4a1a":["D4a1a1"], "D4a1a1":["D4a1a1a"], "D4a1b":["D4a1b1"], "D4a1e":["D4a1e1"], "D4a1f":["D4a1f1"], "D4a2":["D4a2a","D4a2b"], "D4a3":["D4a3a","D4a3b"], "D4a3a":["D4a3a1","D4a3a2"], "D4a3b":["D4a3b1","D4a3b2"], "D4b":["D4b1","D4b2"], "D4b1":["D4b1a","D4b1b'd","D4b1c"], "D4b1a":["D4b1a1","D4b1a2"], "D4b1a1":["D4b1a1a"], "D4b1a2":["D4b1a2a"], "D4b1a2a":["D4b1a2a1","D4b1a2a2"], "D4b1b'd":["D4b1b","D4b1d"], "D4b1b":["D4b1b1","D4b1b2"], "D4b1b1":["D4b1b1a"], "D4b1b1a":["D4b1b1a1"], "D4b1c":["D3"], "D4b2":["D4b2a","D4b2b","D4b2d"], "D4b2a":["D4b2a1","D4b2a2"], "D4b2a2":["D4b2a2a","D4b2a2b"], "D4b2a2a":["D4b2a2a1","D4b2a2a2"], "D4b2b":["D4b2b1","D4b2b2","D4b2b3","D4b2b4","D4b2b5","D4b2b6","D4b2b7"], "D4b2b1":["D4b2b1a","D4b2b1b","D4b2b1c","D4b2b1d"], "D4b2b2":["D4b2b2a","D4b2b2b","D4b2b2c"], "D4b2b2a":["D4b2b2a1"], "D4c":["D4c1","D4c2"], "D4c1":["D4c1a","D4c1b"], "D4c1a":["D4c1a1"], "D4c1b":["D4c1b1","D4c1b2"], "D4c2":["D4c2a","D4c2b","D4c2c"], "D4e":["D4e1'3","D4e2","D4e4","D4e5"], "D4e1'3":["D4e1","D4e3"], "D4e1":["D4e1a","D4e1c","D2"], "D4e1a":["D4e1a1","D4e1a2","D4e1a3"], "D4e1a2":["D4e1a2a"], "D2":["D2a'b","D2c"], "D2a'b":["D2a","D2b"], "D2a":["D2a1","D2a2"], "D2a1":["D2a1a","D2a1b"], "D2b":["D2b1","D2b2"], "D2b1":["D2b1a"], "D4e2":["D4e2a","D4e2b","D4e2c","D4e2d"], "D4e4":["D4e4a","D4e4b"],"D4e4a":["D4e4a1"], "D4e5":["D4e5a","D4e5b"], "D4f":["D4f1"], "D4g":["D4g1","D4g2"], "D4g1":["D4g1a","D4g1b","D4g1c"], "D4g2":["D4g2a","D4g2b"], "D4g2a":["D4g2a1"], "D4g2a1":["D4g2a1a","D4g2a1b","D4g2a1c"], "D4g2b":["D4g2b1"], "D4g2b1":["D4g2b1a"], "D4h":["D4h1","D4h2","D4h3","D4h4"], "D4h1":["D4h1a","D4h1b","D4h1c","D4h1d"], "D4h1a":["D4h1a1","D4h1a2"], "D4h1c":["D4h1c1"], "D4h3":["D4h3a","D4h3b"], "D4h3a":["D4h3a1","D4h3a2","D4h3a3","D4h3a4","D4h3a5","D4h3a6","D4h3a7","D4h3a8","D4h3a9"], "D4h3a1":["D4h3a1a"], "D4h3a1a":["D4h3a1a1","D4h3a1a2"], "D4h3a3":["D4h3a3a"], "D4h4":["D4h4a"], "D4i":["D4i1","D4i2","D4i3"], "D4j":["D4j1","D4j2","D4j3","D4j4","D4j5","D4j6","D4j13","D4j7","D4j8","D4j9","D4j10","D4j12","D4j14","D4j15","D4j16"], "D4j1":["D4j1a","D4j1b"], "D4j1a":["D4j1a1","D4j1a2"], "D4j1a1":["D4j1a1a","D4j1a1b"], "D4j1b":["D4j1b2"], "D4j2":["D4j2a"], "D4j3":["D4j3a"], "D4j3a":["D4j3a1"], "D4j4":["D4j4a"], "D4j5":["D4j5a"], "D4j7":["D4j7a"], "D4o":["D4o1","D4o2"], "D4o1":["D4o1a"], "D4o2":["D4o2a"], "D4o2a":["D4o2a1"], "D4p":["D4p1"], "D4l":["D4l1","D4l2"], "D4l1":["D4l1a"], "D4l1a":["D4l1a1"], "D4l2":["D4l2a","D4l2b"], "D4l2a":["D4l2a1","D4l2a2"], "D4m":["D4m1","D4m2"], "D4m2":["D4m2a"], "D4m2a":["D4m2a1"], "D4m2a1":["D4m2a1a"], "D4n":["D4n1","D4n2"], "D4n1":["D4n1a"], "D4q":["D4q1"], "D4q1":["D4q1a"], "D5":["D5a'b","D5c"], "D5a'b":["D5a","D5b"], "D5a":["D5a1","D5a2","D5a3"], "D5a1":["D5a1a"], "D5a1a":["D5a1a1","D5a1a2"], "D5a2":["D5a2a","D5a2b"], "D5a2a":["D5a2a1","D5a2a2"], "D5a2a1":["D5a2a1a","D5a2a1b"], "D5a2a1a":["D5a2a1a1","D5a2a1a2"], "D5a2a1a1":["D5a2a1a1a"], "D5a2a1b":["D5a2a1b1"], "D5a3":["D5a3a"], "D5a3a":["D5a3a1"], "D5a3a1":["D5a3a1a"], "D5b":["D5b1","D5b2","D5b3","D5b4"], "D5b1":["D5b1a","D5b1b","D5b1c","D5b1d"], "D5b1a":["D5b1a1","D5b1a2"], "D5b1b":["D5b1b1","D5b1b2"], "D5b1c":["D5b1c1"],"D5b1c1":["D5b1c1a"], "D5b3":["D5b3a"], "D5b3a":["D5b3a1"], "D5c":["D5c1","D5c2"], "D5c1":["D5c1a"], "D6":["D6a","D6c"], "D6a":["D6a1","D6a2"], "D6a1":["D6a1a"], "D6c":["D6c1"], "D6c1":["D6c1a"]}

            mt_tree_N = {"N":["R","N1'5","N2","N3","N7","N8","N9","N10","N11","N13","N14","N21","N22","A","O","S","X"], "N1'5":["N1","N5","N5a"],
             "N1":["N1a","N1b"], "N1a":["N1a1'2","N1a3"], "N1a1'2":["N1a1","N1a2"], "N1a1":["N1a1a","N1a1b"], "N1a1a":["N1a1a1","N1a1a2","N1a1a3"], "N1a1a1":["N1a1a1a","N1a1a1b"], "N1a1a1a":["N1a1a1a1","N1a1a1a2","N1a1a1a3"], "N1a1a1a1":["N1a1a1a1a"], "N1a1b":["N1a1b1","I"], "N1a3":["N1a3a"], "N1a3a":["N1a3a1","N1a3a2","N1a3a3"], "N1a3a1":["N1a3a1a"], "N1b":["N1b1","N1b2"], "N1b1":["N1b1a","N1b1b"], "N1b1a":["N1b1a1","N1b1a2","N1b1a3","N1b1a4","N1b1a5","N1b1a6","N1b1a7","N1b1a8"], "N1b1a2":["N1b1a2a","N1b1a2b"], "N1b1a4":["N1b1a4a"], "N1b1a8":["N1b1a8a","N1b1a8b"], "N1b1b":["N1b1b1"],
             "N2":["N2a","W"], "N2a":["N2a1","N2a2"],
             "N3":["N3a","N3b"], "N3a":["N3a1"], "N7":["N7a","N7b"], "N7a":["N7a1","N7a2"],
             "N9":["N9a","N9b","Y"], "N9a":["N9a1'3","N9a2'4'5'11", "N9a6","N9a7","N9a8","N9a9","N9a10"], "N9a1'3":["N9a1","N9a3"], "N9a1":["N9a1a"], "N9a2'4'5'11":["N9a2","N9a4","N9a5","N9a11"], "N9a2":["N9a2a","N9a2c","N9a2d"], "N9a2a":["N9a2a1","N9a2a2","N9a2a3"], "N9a4":["N9a4a","N9a4b"], "N9a4b":["N9a4b1"],"N9a6":["N9a6a","N9a6b"], "N9a10":["N9a10a","N9a10b"], "N9a10a":["N9a10a1","N9a10a2"], "N9a10a2":["N9a10a2a"], "N9b":["N9b1","N9b2","N9b3","N9b4"], "N9b1":["N9b1a","N9b1b","N9b1c"], "N9b1c":["N9b1c1"], "N9b2":["N9b2a"],
             "N10":["N10a","N10b"], "N11":["N11a","N11b"], "N11a":["N11a1","N11a2"], "N21":["N21a"], "N22":["N22a"]}

            mt_tree_A = {"A":["A5","A8","A10","A1","A2","A6","A12","A23","A13","A14","A15","A16","A17","A18","A19","A20","A21","A22","A24","A25","A26","A3","A7","A9","A11"], "A1":["A1a"], "A1a":["A1a1"], "A2":["A2a","A2b","A2c","A2d","A2e","A2ao","A2f","A2g","A2h","A2i","A2j","A2k","A2l","A2m","A2n","A2o","A2ai","A2aj","A2p","A2am","A2q","A2t","A2u","A2v","A2w","A2x","A2y","A2aa","A2ab","A2ac","A2ad","A2ae","A2af","A2ag","A2ah","A2ak","A2al","A2an","A2ap","A2aq","A2r","A2s","A2z"], "A2a":["A2a1","A2a2","A2a3","A2a4","A2a5"], "A2b":["A2b1"], "A2d":["A2d1","A2d2"], "A2d1":["A2d1a"], "A2ao":["A2ao1"], "A2f":["A2f1","A2f2","A2f3"], "A2f1":["A2f1a"], "A2g":["A2g1"], "A2h":["A2h1"], "A2j":["A2j1"], "A2k":["A2k1"], "A2k1":["A2k1a"], "A2p":["A2p1","A2p2"], "A2q":["A2q1"], "A2u":["A2u1","A2u2"], "A2v":["A2v1"], "A2v1":["A2v1a","A2v1b"], "A2w":["A2w1"], "A2ac":["A2ac1"], "A2ad":["A2ad1","A2ad2"], "A2af":["A2af1","A2af2"], "A2af1":["A2af1a","A2af1b"], "A2af1a":["A2af1a1","A2af1a2"], "A2af1b":["A2af1b1","A2af1b2"], "A2af1b1":["A2af1b1a","A2af1b1b"], "A2r":["A2r1"], "A6":["A6a","A6b"], "A12":["A12a"], "A15":["A15a","A15b","A15c"], "A15c":["A15c1"], "A3":["A3a"], "A11":["A11a","A11b"], "A5":["A5a","A5b","A5c"], "A5a":["A5a1","A5a2","A5a3","A5a4","A5a5"], "A5a1":["A5a1a","A5a1b"], "A5a1a":["A5a1a1","A5a1a2"], "A5a1a1":["A5a1a1a","A5a1a1b"], "A5a1a2":["A5a1a2a"], "A5a3":["A5a3a"], "A5b":["A5b1"], "A5b1":["A5b1a","A5b1b","A5b1c"], "A5b1c":["A5b1c1"], "A5c":["A5c1"], "A8":["A8a"], "A8a":["A8a1"]}

            mt_tree_O = {"O":["O1"], "O1":["O1a"]}

            mt_tree_S = {"S":["S1","S2","S3","S4","S5"], "S1":["S1a"]}

            mt_tree_X = {"X":["X1'2'3","X4"], "X1'2'3":["X1'3","X2"], "X1'3":["X1","X3"], "X1":["X1a","X1c"], "X3":["X3a"], "X2":["X2f","X2k","X2p","X2a'j","X2b'd","X2c","X2e","X2g","X2l","X2h","X2i","X2m'n","X2o"], "X2a'j":["X2a","X2j"], "X2a":["X2a1","X2a2"], "X2a1":["X2a1a","X2a1b","X2a1c"], "X2a1a":["X2a1a1"], "X2a1b":["X2a1b1"], "X2a1b1":["X2a1b1a"], "X2b'd":["X2b","X2d"], "X2b":["X2b1","X2b2","X2b3","X2b4","X2b5","X2b6","X2b7","X2b8","X2b9","X2b9","X2b10","X2b11","X2b12","X2b13"], "X2b4":["X2b4a"], "X2b4a":["X2b4a1"], "X2b6":["X2b6a"], "X2b10":["X2b10a"], "X2d":["X2d1","X2d2"], "X2d1":["X2d1a"], "X2c":["X2c1","X2c2"], "X2c1":["X2c1a","X2c1b","X2c1c","X2c1d","X2c1e"], "X2c1c":["X2c1c1"], "X2e":["X2e1","X2e2"], "X2e1":["X2e1a","X2e1b"], "X2e1a":["X2e1a1"], "X2e2":["X2e2a","X2e2b","X2e2c"], "X2e2a":["X2e2a1","X2e2a2"], "X2e2b":["X2e2b1"], "X2e2c":["X2e2c1"], "X2i":["X2i1"], "X2m'n":["X2m","X2n"], "X2m":["X2m1","X2m2"], "X2o":["X2o1"], "X2f":["X2f1"], "X2p":["X2p1"]}

            mt_tree_I = {"I":["I1","I2'3","I4","I5","I6","I7"], "I1":["I1a","I1b","I1c","I1d","I1e","I1f"], "I1a":["I1a1"], "I1a1":["I1a1a","I1a1b","I1a1c","I1a1d","I1a1e"], "I1a1a":["I1a1a1","I1a1a2","I1a1a3"], "I1a1a3":["I1a1a3a"], "I1c":["I1c1"], "I1c1":["I1c1a"], "I2'3":["I2","I3"], "I2":["I2a","I2b","I2c","I2d","I2e","I2f"], "I2a":["I2a1","I2a2","I2a3"], "I2a1":["I2a1a"], "I3":["I3a","I3b","I3c","I3d"], "I3a":["I3a1"], "I3d":["I3d1"], "I4":["I4a","I4b"], "I4a":["I4a1","I4a2"], "I5":["I5a","I5b","I5c"], "I5a":["I5a1","I5a2","I5a3","I5a4"], "I5a1":["I5a1a","I5a1b","I5a1c"], "I5a2":["I5a2a"], "I5b":["I5b1"], "I5c":["I5c1"], "I6":["I6a","I6b"]}

            mt_tree_R = {"R":["R0","R1","R2'JT","R5","R6","R7","R8","R9","R11'B6","B4'5","R24","R12'21","R14","R22","R23","R30","R31","R32","P","U"],
             "R0":["R0a'b","HV"], "R0a'b":["R0a","R0b"], "R0a":["R0a1","R0a2'3","R0a4"], "R0a1":["R0a1a","R0a1b"], "R0a1a":["R0a1a1","R0a1a2","R0a1a3","R0a1a4"], "R0a1a1":["R0a1a1a"], "R0a2'3":["R0a2","R0a3"], "R0a2":["R0a2a","R0a2b","R0a2c","R0a2d","R0a2e","R0a2f","R0a2g","R0a2h","R0a2i","R0a2j","R0a2k","R0a2l","R0a2m","R0a2n"], "R0a2a":["R0a2a1"], "R0a2f":["R0a2f1"], "R0a2f1":["R0a2f1a","R0a2f1b"], "R0a2k":["R0a2k1"], "R0a3":["R0a3a"], 
             "R1":["R1a","R1b"], "R1a":["R1a1"], "R1a1":["R1a1a","R1a1b","R1a1c"], "R1a1a":["R1a1a1","R1a1a2"], "R1a1a1":["R1a1a1a"], "R1b":["R1b1"],
             "R2'JT":["R2","JT"], "R2":["R2a","R2b","R2c","R2d"], "R2b":["R2b1"],
             "R5":["R5a"], "R5a":["R5a1","R5a2"], "R5a1":["R5a1a"], "R5a2":["R5a2a","R5a2b"], "R5a2b":["R5a2b1","R5a2b2","R5a2b3","R5a2b4"], "R6":["R6a","R6b"], "R6a":["R6a1","R6a2"], "R7":["R7a'b"], "R7a'b":["R7a","R7b"], "R7a":["R7a1"], "R7a1":["R7a1a","R7a1b"], "R7a1b":["R7a1b1","R7a1b2"], "R7b":["R7b1","R7b2"], "R7b1":["R7b1a"], "R7b1a":["R7b1a1"], "R8":["R8a","R8b"], "R8a":["R8a1","R8a2"], "R8a1":["R8a1a","R8a1b"], "R8a1a":["R8a1a1","R8a1a2","R8a1a3"], "R8a1a1":["R8a1a1a","R8a1a1b","R8a1a1c","R8a1a1d"], "R8a1a1a":["R8a1a1a1","R8a1a1a2"], "R8a1a1a1":["R8a1a1a1a"], "R8a1a2":["R8a1a2a"], "R8b":["R8b1","R8b2"], "R8b1":["R8b1a"],
             "R9":["R9b","R9c","F"], "R9b":["R9b1","R9b2"], "R9b1":["R9b1a","R9b1b"], "R9b1a":["R9b1a1","R9b1a2","R9b1a3"], "R9b1a1":["R9b1a1a"], "R9b1a2":["R9b1a2a","R9b1a2b"], "R9c":["R9c1"], "R9c1":["R9c1a","R9c1b"], "R9c1a":["R9c1a1","R9c1a2","R9c1a3"], "R9c1b":["R9c1b1","R9c1b2"],
             "R11'B6":["R11","B6"], "R11":["R11a","R11b"], "R11b":["R11b1"], "R11b1":["R11b1a","R11b1b"], "R24":["R24a"], "R12'21":["R12","R21"], "R30":["R30a","R30b"], "R30a":["R30a1"], "R30a1":["R30a1a","R30a1b","R30a1c"], "R30a1b":["R30a1b1"], "R30b":["R30b1","R30b2"], "R30b2":["R30b2a"], "R31":["R31a","R31b"], "R31a":["R31a1"]}

            mt_tree_W = {"W":["W1","W3","W4","W5","W6","W7","W8","W9"], "W1":["W1a","W1b","W1c","W1i","W1d","W1e","W1f","W1g","W1h"], "W1b":["W1b1"], "W1c":["W1c1"], "W1e":["W1e1"], "W1e1":["W1e1a"], "W1h":["W1h1"],"W3":["W3a","W3b"], "W3a":["W3a1","W3a2"], "W3a1":["W3a1a","W3a1b","W3a1c","W3a1d"], "W3a1a":["W3a1a1","W3a1a2","W3a1a3"], "W3b":["W3b1"], "W4":["W4a","W4b","W4c","W4d"], "W4a":["W4a1"], "W5":["W5a","W5b"], "W5a":["W5a1","W5a2"], "W5a1":["W5a1a"], "W5a1a":["W5a1a1"], "W5a1a1":["W5a1a1a"], "W5a2":["W5a2b"], "W5b":["W5b1"], "W5b1":["W5b1a"], "W6":["W6a","W6b","W6c","W6d"], "W6b":["W6b1"], "W6c":["W6c1"], "W6c1":["W6c1a"]}

            mt_tree_Y = {"Y":["Y1","Y2"]}

            mt_tree_P = {"P":["P1","P2'10","P8","P3","P4","P5","P6","P7","P9"], "P1":["P1d","P1e","P1f"], "P1d":["P1d1","P1d2"], "P1d1":["P1d1a"], "P1d2":["P1d2a"], "P2'10":["P2","P10"], "P3":["P3a","P3b"], "P3b":["P3b1"], "P4":["P4a","P4b"], "P4a":["P4a1"], "P4b":["P4b1"], "P9":["P9a"]}

            mt_tree_U = {"U":["U1","U5","U6","U2'3'4'7'8'9"],
             "U1":["U1a","U1b"], "U1a":["U1a1","U1a2","U1a3"], "U1a1":["U1a1a","U1a1b","U1a1c","U1a1d"], "U1a1a":["U1a1a1","U1a1a2","U1a1a3"], "U1a1a1":["U1a1a1a"], "U1a1c":["U1a1c1"], "U1a1c1":["U1a1c1a","U1a1c1b","U1a1c1c","U1a1c1d"], "U1a1c1c":["U1a1c1c1"], "U1a1c1d":["U1a1c1d1"], "U1b":["U1b1","U1b2","U1b3"],
             "U5":["U5a'b"], "U5a'b":["U5a","U5b"], "U5a":["U5a1","U5a2"], "U5a1":["U5a1a","U5a1g","U5a1b","U5a1c","U5a1d","U5a1e","U5a1f","U5a1h","U5a1i","U5a1j"], "U5a1a":["U5a1a1","U5a1a2"], "U5a1a1":["U5a1a1a","U5a1a1b","U5a1a1h","U5a1a1c","U5a1a1d","U5a1a1e","U5a1a1g","U5a1a1i"], "U5a1a1d":["U5a1a1d1"], "U5a1a2":["U5a1a2a","U5a1a2b"], "U5a1a2a":["U5a1a2a1"], "U5a1a2a1":["U5a1a2a1a"], "U5a1a2b":["U5a1a2b1"], "U5a1g":["U5a1g1","U5a1g2"], "U5a1b":["U5a1b1","U5a1b2", "U5a1b3","U5a1b4"], "U5a1b1":["U5a1b1a","U5a1b1b","U5a1b1c","U5a1b1d","U5a1b1e","U5a1b1f","U5a1b1g","U5a1b1h"], "U5a1b1a":["U5a1b1a1","U5a1b1a2"], "U5a1b1b":["U5a1b1b1"], "U5a1b1c":["U5a1b1c1","U5a1b1c2"], "U5a1b1d":["U5a1b1d1"], "U5a1b3":["U5a1b3a"], "U5a1b3a":["U5a1b3a1"], "U5a1c":["U5a1c1","U5a1c2"], "U5a1c1":["U5a1c1a"], "U5a1c2":["U5a1c2a"], "U5a1c2a":["U5a1c2a1"], "U5a1d":["U5a1d1","U5a1d2"], "U5a1d2":["U5a1d2a","U5a1d2b"], "U5a1d2a":["U5a1d2a1"], "U5a1f":["U5a1f1","U5a1f2"], "U5a1f1":["U5a1f1a"], "U5a1f1a":["U5a1f1a1"], "U5a1i":["U5a1i1"], "U5a2":["U5a2a","U5a2b","U5a2c","U5a2d","U5a2e"], "U5a2a":["U5a2a1","U5a2a2"], "U5a2a1":["U5a2a1a","U5a2a1b","U5a2a1c","U5a2a1d","U5a2a1e"], "U5a2a1b":["U5a2a1b1"], "U5a2a2":["U5a2a2a"], "U5a2b":["U5a2b1","U5a2b2","U5a2b3","U5a2b4","U5a2b5"], "U5a2b1":["U5a2b1a","U5a2b1b","U5a2b1c","U5a2b1d"], "U5a2b2":["U5a2b2a"], "U5a2b2a":["U5a2b2a1"], "U5a2b3":["U5a2b3a","U5a2b3a1"], "U5a2b4":["U5a2b4a"], "U5a2c":["U5a2c1","U5a2c2","U5a2c3","U5a2c4"], "U5a2c3":["U5a2c3a"], "U5a2d":["U5a2d1"], "U5a2d1":["U5a2d1a"], "U5b":["U5b1","U5b2","U5b3"], "U5b1":["U5b1a","U5b1b","U5b1c","U5b1e","U5b1h","U5b1d","U5b1f","U5b1g","U5b1i"], "U5b1b":["U5b1b1","U5b1b2"], "U5b1b1":["U5b1b1a","U5b1b1d","U5b1b1f","U5b1b1b","U5b1b1e","U5b1b1g"], "U5b1b1a":["U5b1b1a1","U5b1b1a2","U5b1b1a3"], "U5b1b1a1":["U5b1b1a1a","U5b1b1a1b"], "U5b1b1a1a":["U5b1b1a1a1"], "U5b1b1g":["U5b1b1g1"], "U5b1b1g1":["U5b1b1g1a"], "U5b1b2":["U5b1b2a","U5b1b2b"], "U5b1c":["U5b1c1","U5b1c2"], "U5b1c1":["U5b1c1a"], "U5b1c1a":["U5b1c1a1"], "U5b1c2":["U5b1c2a","U5b1c2b"], "U5b1e":["U5b1e1"], "U5b1e1":["U5b1e1a"], "U5b1d":["U5b1d1","U5b1d2"], "U5b1d1":["U5b1d1a","U5b1d1b","U5b1d1c"], "U5b1f":["U5b1f1"], "U5b1f1":["U5b1f1a"], "U5b2":["U5b2a","U5b2b","U5b2c"], "U5b2a":["U5b2a1","U5b2a2","U5b2a3","U5b2a4","U5b2a5","U5b2a6"], "U5b2a1":["U5b2a1a","U5b2a1b"], "U5b2a1a":["U5b2a1a1a","U5b2a1a1b","U5b2a1a1d"], "U5b2a2":["U5b2a2a","U5b2a2b","U5b2a2c"], "U5b2a2a":["U5b2a2a1","U5b2a2a2"], "U5b2a2b":["U5b2a2b1"], "U5b2a3":["U5b2a3a"], "U5b2a4":["U5b2a4a"], "U5b2a5":["U5b2a5a"], "U5b2b":["U5b2b1","U5b2b2","U5b2b3","U5b2b4","U5b2b5"], "U5b2b1":["U5b2b1a","U5b2b1b"], "U5b2b1a":["U5b2b1a1","U5b2b1a2"], "U5b2b3":["U5b2b3a","U5b2b3b"], "U5b2b3a":["U5b2b3a1"], "U5b2b3a1":["U5b2b3a1a"], "U5b2b4":["U5b2b4a"], "U5b2c":["U5b2c1","U5b2c2"], "U5b2c2":["U5b2c2a","U5b2c2b"], "U5b3":["U5b3a","U5b3b","U5b3c","U5b3d","U5b3e","U5b3f","U5b3g","U5b3h"], "U5b3a":["U5b3a1","U5b3a2"], "U5b3a1":["U5b3a1a","U5b3a1b"], "U5b3b":["U5b3b1","U5b3b2"],
             "U6":["U6a'b'd","U6c"],"U6a'b'd":["U6a","U6b","U6d"], "U6a":["U6a1","U6a2","U6a8","U6a3","U6a4","U6a5","U6a6","U6a7"], "U6a1":["U6a1a","U6a1b"], "U6a1a":["U6a1a1","U6a1a2"], "U6a1b":["U6a1b1","U6a1b2","U6a1b3","U6a1b4"], "U6a1b1":["U6a1b1a","U6a1b1b"], "U6a2":["U6a2a","U6a2b"], "U6a2a":["U6a2a1","U6a2a2"], "U6a2a2":["U6a2a2a"], "U6a2b":["U6a2b1"], "U6a8":["U6a8a","U6a8b"], "U6a3":["U6a3a","U6a3b","U6a3e","U6a3f","U6a3c","U6a3d"], "U6a3a":["U6a3a1","U6a3a2"], "U6a3a1":["U6a3a1a"], "U6a3a2":["U6a3a2a"], "U6a3b":["U6a3b1"], "U6a3f":["U6a3f1","U6a3f2"],"U6a3d":["U6a3d1"], "U6a3d1":["U6a3d1a"], "U6a5":["U6a5a","U6a5b","U6a5c"], "U6a5a":["U6a5a1"], "U6a6":["U6a6a","U6a6b"], "U6a6a":["U6a6a1"], "U6a6b":["U6a6b1","U6a6b2"], "U6a7":["U6a7a","U6a7b","U6a7c"],"U6a7a":["U6a7a1","U6a7a2"], "U6a7a1":["U6a7a1a","U6a7a1b","U6a7a1c"], "U6a7a2":["U6a7a2a"], "U6a7b":["U6a7b1"], "U6a7c":["U6a7c1"], "U6b":["U6b1","U6b2","U6b3"], "U6b1":["U6b1a","U6b1b"], "U6b1a":["U6b1a1"], "U6b3":["U6b3a"], "U6d":["U6d1","U6d2","U6d3"], "U6d1":["U6d1a","U6d1b"], "U6d3":["U6d3a"], "U6c":["U6c1","U6c2"],
             "U2'3'4'7'8'9":["U2","U3","U4","U4'9","U7","U8"], 
             "U2":["U2a","U2b", "U2c'd", "U2e"], "U2a":["U2a1","U2a2"], "U2a1":["U2a1a","U2a1b"], "U2b":["U2b1","U2b2"], "U2b1":["U2b1a"], "U2c'd":["U2c","U2d"], "U2c":["U2c1"], "U2c1":["U2c1a","U2c1b"], "U2d":["U2d1","U2d2","U2d3"], "U2d2":["U2d2a"], "U2e":["U2e1'2'3"], "U2e1'2'3":["U2e1","U2e2","U2e3"], "U2e1":["U2e1a","U2e1b","U2e1c","U2e1d","U2e1e","U2e1f","U2e1g","U2e1h"], "U2e1a":["U2e1a1"], "U2e1a1":["U2e1a1a","U2e1a1b","U2e1a1c"], "U2e1b":["U2e1b1","U2e1b2"], "U2e1c":["U2e1c1"], "U2e1f":["U2e1f1"], "U2e2":["U2e2a"], "U2e2a":["U2e2a1"], "U2e2a1":["U2e2a1a","U2e2a1b","U2e2a1c","U2e2a1d"], "U2e2a1a":["U2e2a1a1","U2e2a1a2"], "U2e3":["U2e3a"],
             "U3":["U3a'c","U3b"], "U3a'c":["U3a","U3c"], "U3a":["U3a1","U3a2","U3a3"], "U3a1":["U3a1a","U3a1b","U3a1c"], "U3a1a":["U3a1a1"], "U3a1c":["U3a1c1"], "U3a2":["U3a2a"], "U3a2a":["U3a2a1"], "U3a2a1":["U3a2a1a"], "U3b":["U3b1","U3b2","U3b3"], "U3b1":["U3b1a","U3b1b"], "U3b1a":["U3b1a1"], "U3b2":["U3b2a","U3b2b","U3b2c"], "U3b2a":["U3b2a1"],"U3b2a1":["U3b2a1a"], "U4'9":["U4","U9"], "U4":["U4a","U4b","U4c","U4d"], "U4a":["U4a1","U4a2","U4a3"], "U4a1":["U4a1a","U4a1b","U4a1c","U4a1d","U4a1e"], "U4a1a":["U4a1a1","U4a1a2","U4a1a3"], "U4a1b":["U4a1b1","U4a1b2"], "U4a1b1":["U4a1b1a"], "U4a2":["U4a2a","U4a2b","U4a2c","U4a2d","U4a2e","U4a2f","U4a2g","U4a2h"], "U4a2a":["U4a2a1","U4a2a2","U4a2a3"], "U4a2c":["U4a2c1"], "U4a2h":["U4a2h1"], "U4a3":["U4a3a"], "U4b":["U4b1","U4b2","U4b3"], "U4b1":["U4b1a","U4b1b"], "U4b1a":["U4b1a1","U4b1a2","U4b1a3","U4b1a4"], "U4b1a1":["U4b1a1a"], "U4b1a1a":["U4b1a1a1"], "U4b1a2":["U4b1a2a","U4b1a2b"], "U4b1a3":["U4b1a3a"], "U4b1b":["U4b1b1","U4b1b2"], "U4b1b1":["U4b1b1a","U4b1b1b","U4b1b1c","U4b1b1d"], "U4b2":["U4b2a"], "U4b2a":["U4b2a1"], "U4b2a1":["U4b2a1a"], "U4c":["U4c1","U4c2"], "U4c1":["U4c1a"], "U4c2":["U4c2a"], "U4d":["U4d1","U4d2","U4d3"], "U4d1":["U4d1a","U4d1b"], "U4d1a":["U4d1a1"], "U4d1a1":["U4d1a1a"], "U9":["U9a","U9b"], "U9a":["U9a1"], "U9b":["U9b1"],
             "U7":["U7a","U7b"], "U7a":["U7a1","U7a2","U7a3","U7a4","U7a5"], "U7a1":["U7a1a"], "U7a2":["U7a2a"], "U7a3":["U7a3a","U7a3b"], "U7a4":["U7a4a"], "U7a4a":["U7a4a1"], "U7a4a1":["U7a4a1a"], "U7b":["U7b1","U7b2"],
             "U8":["U8a","U8b'c"], "U8a":["U8a1","U8a2"], "U8a1":["U8a1a","U8a1b"], "U8a1a":["U8a1a1","U8a1a2","U8a1a3","U8a1a4"], "U8a1a1":["U8a1a1a","U8a1a1b"], "U8a1a1a":["U8a1a1a1"], "U8a1a1b":["U8a1a1b1"], "U8b'c":["U8b","U8c"], "U8b":["U8b1","K"], "U8b1":["U8b1a","U8b1b"], "U8b1a":["U8b1a1","U8b1a2"], "U8b1a2":["U8b1a2a","U8b1a2b"], "U8b1b":["U8b1b1","U8b1b2"]}

            mt_tree_J = {"J":["J1","J2"], 
             "J1":["J1b","J1c","J1d"], "J1b":["J1b1","J1b2","J1b3","J1b4","J1b5","J1b6","J1b7","J1b8","J1b9"], "J1b1":["J1b1a","J1b1b"], "J1b1a":["J1b1a1","J1b1a2","J1b1a3"], "J1b1a1":["J1b1a1a","J1b1a1b","J1b1a1c","J1b1a1d","J1b1a1e"], "J1b1a2":["J1b1a2a","J1b1a2b"], "J1b1b":["J1b1b1","J1b1b2","J1b1b3"], "J1b1b1":["J1b1b1a","J1b1b1b","J1b1b1c"], "J1b2":["J1b2a"], "J1b3":["J1b3a","J1b3b"], "J1b3b":["J1b3b1"], "J1b4":["J1b4a"], "J1b4a":["J1b4a1","J1b4a2"], "J1b5":["J1b5a"], "J1b5a":["J1b5a1"], "J1b6":["J1b6a","J1b6b"], "J1b7":["J1b7a"], "J1c":["J1c1","J1c2","J1c3","J1c4","J1c5","J1c6","J1c7","J1c12","J1c13","J1c14","J1c8","J1c9","J1c10","J1c11","J1c15","J1c16","J1c17"], "J1c1":["J1c1a","J1c1b","J1c1c","J1c1d","J1c1e","J1c1f","J1c1g","J1c1h"], "J1c1b":["J1c1b1","J1c1b2"], "J1c1b1":["J1c1b1a"], "J1c1b1a":["J1c1b1a1"], "J1c1b2":["J1c1b2a"], "J1c1g":["J1c1g1"], "J1c2":["J1c2a","J1c2b","J1c2c","J1c2d","J1c2e","J1c2f","J1c2g","J1c2h","J1c2i","J1c2j","J1c2k","J1c2l","J1c2m","J1c2n","J1c2o","J1c2p","J1c2q","J1c2r","J1c2s","J1c2t"], "J1c2a":["J1c2a1","J1c2a2","J1c2a3"], "J1c2a1":["J1c2a1a"], "J1c2b":["J1c2b1","J1c2b2","J1c2b3","J1c2b4","J1c2b5"], "J1c2c":["J1c2c1","J1c2c2","J1c2c3"], "J1c2c1":["J1c2c1a"], "J1c2c2":["J1c2c2a"], "J1c2e":["J1c2e1","J1c2e2"], "J1c2m":["J1c2m1"], "J1c2n":["J1c2n1"], "J1c2q":["J1c2q1"], "J1c2s":["J1c2s1"], "J1c3":["J1c3a","J1c3b","J1c3c","J1c3d","J1c3e","J1c3f","J1c3g","J1c3h","J1c3i","J1c3j","J1c3k","J1c3m","J1c3n'o'p"], "J1c3a":["J1c3a1","J1c3a2"], "J1c3b":["J1c3b1","J1c3b2"], "J1c3b1":["J1c3b1a"], "J1c3c":["J1c3c1","J1c3c2"], "J1c3e":["J1c3e1","J1c3e2"], "J1c4":["J1c4b","J1c4c"], "J1c5":["J1c5a","J1c5b","J1c5c","J1c5d","J1c5e","J1c5f"], "J1c5a":["J1c5a1"], "J1c5c":["J1c5c1"], "J1c6":["J1c6a"], "J1c7":["J1c7a"], "J1c12":["J1c12a","J1c12b"], "J1c8":["J1c8a","J1c8b"], "J1c8a":["J1c8a1"], "J1c8a1":["J1c8a1a"], "J1c10":["J1c10a"], "J1c11":["J1c11a"], "J1c15":["J1c15a","J1c15b"], "J1c15a":["J1c15a1"], "J1c17":["J1c17a"], "J1d":["J1d1","J1d2","J1d3","J1d4","J1d5","J1d6"], "J1d1":["J1d1a","J1d1b"], "J1d1a":["J1d1a1"], "J1d1a1":["J1d1a1a"], "J1d1b":["J1d1b1"], "J1d2":["J1d2a"], "J1d3":["J1d3a"], "J1d3a":["J1d3a1","J1d3a2"], "J1d5":["J1d5a"], "J1d6":["J1d6a"], 
             "J2":["J2a","J2b"], "J2a":["J2a1","J2a2"], "J2a1":["J2a1a"], "J2a1a":["J2a1a1","J2a1a2"], "J2a1a1":["J2a1a1a","J2a1a1b","J2a1a1c","J2a1a1d","J2a1a1e"], "J2a1a1a":["J2a1a1a1","J2a1a1a2","J2a1a1a3"], "J2a1a1a2":["J2a1a1a2a"], "J2a1a2":["J2a1a2a"], "J2a1a2a":["J2a1a2a1"], "J2a1a2a1":["J2a1a2a1a"], "J2a2":["J2a2a","J2a2b","J2a2c","J2a2d","J2a2e"], "J2a2a":["J2a2a1","J2a2a2"], "J2a2a1":["J2a2a1a"], "J2a2a1a":["J2a2a1a1"], "J2a2b":["J2a2b1","J2a2b2","J2a2b3"], "J2a2b1":["J2a2b1a"], "J2a2c":["J2a2c1"], "J2b":["J2b1","J2b2"], "J2b1":["J2b1a","J2b1b","J2b1c","J2b1d","J2b1e", "J2b1f","J2b1g","J2b1h"], "J2b1a":["J2b1a1","J2b1a2","J2b1a3","J2b1a4","J2b1a5","J2b1a6"], "J2b1a1":["J2b1a1a"], "J2b1a2":["J2b1a2a"], "J2b1b":["J2b1b1"], "J2b1c":["J2b1c1"], "J2b1e":["J2b1e1"]}

            mt_tree_T = {"T":["T1","T2","T3"],
             "T1":["T1a","T1b"], "T1a":["T1a1'3","T1a5","T1a6","T1a7","T1a8","T1a9","T1a10", "T1a2","T1a4","T1a11","T1a12","T1a13"], "T1a1'3":["T1a1","T1a3"], "T1a1":["T1a1a","T1a1b","T1a1c","T1a1d","T1a1e","T1a1f","T1a1g","T1a1h","T1a1i","T1a1j","T1a1k","T1a1l","T1a1m","T1a1n","T1a1p","T1a1q","T1a1r"], "T1a1a":["T1a1a1"], "T1a1b":["T1a1b1"], "T1a1k":["T1a1k1","T1a1k2"], "T1a1m":["T1a1m1"], "T1a3":["T1a3a"], "T1a2":["T1a2a","T1a2b"], "T1a5":["T1a5a"], "T1a8":["T1a8a"], "T1a10":["T1a10a"], "T1b":["T1b1","T1b2","T1b3","T1b4"],
             "T2":["T2a","T2b","T2c","T2d","T2e","T2m","T2f","T2g","T2h","T2i","T2j","T2k","T2l","T2n"], "T2a":["T2a1","T2a2","T2a3"], "T2a1":["T2a1a","T2a1b"], "T2a1a":["T2a1a1","T2a1a2","T2a1a3","T2a1a5","T2a1a6","T2a1a7","T2a1a8"], "T2a1a3":["T2a1a3a"], "T2a1b":["T2a1b1","T2a1b2"], "T2a1b1":["T2a1b1a"], "T2a1b1a":["T2a1b1a1","T2a1b1a2"], "T2a1b1a1":["T2a1b1a1a","T2a1b1a1b"], "T2a1b1a1a":["T2a1b1a1a1","T2a1b1a1a2"], "T2a1b1a1b":["T2a1b1a1b1"], "T2a1b2":["T2a1b2a","T2a1b2b"], "T2a2":["T2a2a"], "T2b":["T2b1","T2b2","T2b3","T2b4","T2b5","T2b6","T2b7","T2b8","T2b9","T2b11","T2b13","T2b15","T2b16","T2b17","T2b19","T2b21","T2b22","T2b23","T2b24","T2b25","T2b26","T2b27","T2b28","T2b29","T2b30","T2b31","T2b32","T2b33","T2b34","T2b35","T2b36","T2b37"], "T2b2":["T2b2b"], "T2b2b":["T2b2b1"], "T2b3":["T2b3a","T2b3b","T2b3c","T2b3d","T2b3e"], "T2b3a":["T2b3a1"], "T2b4":["T2b4a","T2b4b","T2b4c","T2b4d","T2b4e","T2b4f","T2b4g","T2b4h"], "T2b4a":["T2b4a1"], "T2b5":["T2b5a"], "T2b5a":["T2b5a1"], "T2b6":["T2b6a","T2b6b"], "T2b7":["T2b7a"], "T2b7a":["T2b7a1","T2b7a2","T2b7a3"], "T2b13":["T2b13a","T2b13b"], "T2b17":["T2b17a"], "T2b19":["T2b19b"], "T2b21":["T2b21a","T2b21b"], "T2b23":["T2b23a"], "T2b24":["T2b24a"], "T2c":["T2c1"], "T2c1":["T2c1a","T2c1c","T2c1d","T2c1e","T2c1f"], "T2c1a":["T2c1a1","T2c1a2","T2c1a3"], "T2c1c":["T2c1c1","T2c1c2"], "T2c1d":["T2c1d1","T2c1d2"], "T2c1d1":["T2c1d1a"], "T2c1d2":["T2c1d2a"], "T2d":["T2d1"], "T2d1":["T2d1a","T2d1b"], "T2d1b":["T2d1b1","T2d1b2"], "T2e":["T2e1","T2e2","T2e5","T2e6","T2e7"], "T2e1":["T2e1a","T2e1b"], "T2e1a":["T2e1a1"], "T2e1a1":["T2e1a1a","T2e1a1b"], "T2e1a1b":["T2e1a1b1"], "T2f":["T2f1","T2f2","T2f3","T2f4","T2f5","T2f6","T2f7","T2f8"], "T2f1":["T2f1a"], "T2f1a":["T2f1a1"], "T2f7":["T2f7a"], "T2f8":["T2f8a"], "T2g":["T2g1","T2g2"], "T2g1":["T2g1a","T2g1b"], "T2g1a":["T2g1a1"], "T2g2":["T2g2a"], "T2h":["T2h1","T2h2"], "T2i":["T2i1","T2i2"], "T2j":["T2j1"]}

            mt_tree_B = {"B6":["B6a"], "B6a":["B6a1"], "B6a1":["B6a1a"],
             "B4'5":["B4","B5"],
             "B4":["B4a","B4b'd'e'j","B4c","B4f","B4g","B4h","B4i","B4k","B4m"], "B4a":["B4a1","B4a2","B4a3","B4a4","B4a5"], "B4a1":["B4a1a","B4a1b","B4a1e","B4a1c","B4a1d"], "B4a1a":["B4a1a1","B4a1a2","B4a1a3","B4a1a4","B4a1a5","B4a1a6","B4a1a7"], "B4a1a1":["B4a1a1a","B4a1a1b","B4a1a1c","B4a1a1d","B4a1a1e","B4a1a1f","B4a1a1g","B4a1a1h","B4a1a1i","B4a1a1j","B4a1a1k","B4a1a1n","B4a1a1m","B4a1a1o", "B4a1a1p","B4a1a1q","B4a1a1r","B4a1a1s","B4a1a1t","B4a1a1u","B4a1a1v","B4a1a1w","B4a1a1x","B4a1a1y","B4a1a1z","B4a1a1aa","B4a1a1ab","B4a1a1ac","B4a1a1ad","B4a1a1ae","B4a1a1af"], "B4a1a1a":["B4a1a1a1","B4a1a1a2","B4a1a1a3","B4a1a1a4","B4a1a1a5","B4a1a1a6","B4a1a1a7","B4a1a1a8","B4a1a1a9","B4a1a1a10","B4a1a1a11","B4a1a1a12","B4a1a1a13","B4a1a1a14","B4a1a1a15","B4a1a1a16","B4a1a1a17","B4a1a1a18","B4a1a1a19","B4a1a1a20","B4a1a1a21","B4a1a1a22","B4a1a1a23"], "B4a1a1a1":["B4a1a1a1a","B4a1a1a1b","B4a1a1a1c","B4a1a1a1d"], "B4a1a1a1a":["B4a1a1a1a1"], "B4a1a1a2":["B4a1a1a2a","B4a1a1a2b"], "B4a1a1a11":["B4a1a1a11a","B4a1a1a11b"], "B4a1a1k":["B4a1a1k1"], "B4a1a1m":["B4a1a1m1"], "B4a1a3":["B4a1a3a"], "B4a1a3a":["B4a1a3a1"], "B4a1a3a1":["B4a1a3a1a"], "B4a1a5":["B4a1a5a"], "B4a1a6":["B4a1a6a"], "B4a1b":["B4a1b1"], "B4a1b1":["B4a1b1a"], "B4a1c":["B4a1c1","B4a1c2","B4a1c3","B4a1c4","B4a1c5"], "B4a1c1":["B4a1c1a"], "B4a1c1a":["B4a1c1a1"], "B4a1c3":["B4a1c3a","B4a1c3b"], "B4a2":["B4a2a","B4a2b"], "B4a2a":["B4a2a1","B4a2a2","B4a23"], "B4a2b":["B4a2b1"], "B4a2b1":["B4a2b1a"], "B4g":["B4g1","B4g2"], "B4g1":["B4g1a","B4g1b"], "B4h":["B4h1"], "B4i":["B4i1"], "B4b'd'e'j":["B4b","B4d","B4e","B4j"], "B4b":["B2","B4b1"], "B2":["B2a","B2b","B2c","B2d","B2e","B2f","B2g","B2h","B2i","B2j","B2k","B2l","B2m","B2n","B2o","B2p","B2q","B2r","B2s","B2t","B2u","B2v","B2w","B2x","B2y"], "B2a":["B2a1","B2a2","B2a3","B2a4","B2a5"], "B2a1":["B2a1a","B2a1b"], "B2a1a":["B2a1a1"], "B2a4":["B2a4a"], "B2a4a":["B2a4a1"], "B2b":["B2b1","B2b2","B2b3","B2b4"], "B2b2":["B2b2a"], "B2b3":["B2b3a"], "B2c":["B2c1","B2c2"], "B2c1":["B2c1a","B2c1b","B2c1c"], "B2c2":["B2c2a","B2c2b"], "B2g":["B2g1","B2g2"], "B2i":["B2i1","B2i2"], "B2i2":["B2i2a","B2i2b"], "B2i2a":["B2i2a1"], "B2i2a1":["B2i2a1a","B2i2a1b"], "B2i2b":["B2i2b1"], "B2o":["B2o1"], "B2o1":["B2o1a"], "B2y":["B2y1"], "B4b1":["B4b1a","B4b1b'c"], "B4b1a":["B4b1a1","B4b1a2","B4b1a3"], "B4b1a1":["B4b1a1a","B4b1a1b","B4b1a1c"], "B4b1a2":["B4b1a2a","B4b1a2b","B4b1a2c","B4b1a2d","B4b1a2e","B4b1a2f","B4b1a2g","B4b1a2h","B4b1a2i"], "B4b1a2b":["B4b1a2b1","B4b1a2b2"], "B4b1a2g":["B4b1a2g1"], "B4b1a3":["B4b1a3a"], "B4b1b'c":["B4b1b","B4b1c"], "B4b1c":["B4b1c1","B4b1c2"], "B4d":["B4d1'2'3","B4d4"], "B4d1'2'3":["B4d1","B4d2","B4d3"], "B4d1":["B4d1a"], "B4d3":["B4d3a"], "B4d3a":["B4d3a1"], "B4c":["B4c1","B4c2"], "B4c1":["B4c1a'b","B4c1c"], "B4c1a'b":["B4c1a","B4c1b"], "B4c1a":["B4c1a1","B4c1a2"], "B4c1a1":["B4c1a1a","B4c1a1b","B4c1a1c"], "B4c1a1a":["B4c1a1a1","B4c1a1a2"], "B4c1a1a1":["B4c1a1a1a"], "B4c1a2":["B4c1a2a"], "B4c1b":["B4c1b1","B4c1b2"], "B4c1b1":["B4c1b1a"], "B4c1b2":["B4c1b2a","B4c1b2b","B4c1b2c"], "B4c1b2a":["B4c1b2a1","B4c1b2a2"], "B4c1b2a2":["B4c1b2a2a","B4c1b2a2b"], "B4c1b2c":["B4c1b2c1","B4c1b2c2"], "B4c1c":["B4c1c1"], "B4c2":["B4c2a","B4c2b","B4c2c"], "B4f":["B4f1"], 
             "B5":["B5a","B5b"], "B5a":["B5a1","B5a2"], "B5a1":["B5a1a","B5a1b","B5a1c","B5a1d"], "B5a1a":["B5a1a1"], "B5a1b":["B5a1b1"], "B5a1c":["B5a1c1","B5a1c2"], "B5a1c1":["B5a1c1a"], "B5a1c1a":["B5a1c1a1"], "B5a2":["B5a2a"], "B5a2a":["B5a2a1","B5a2a2"], "B5a2a1":["B5a2a1a","B5a2a1b"], "B5a2a2":["B5a2a2a","B5a2a2b"], "B5a2a2a":["B5a2a2a1","B5a2a2a2"], "B5a2a2b":["B5a2a2b1","B5a2a2b2"], "B5a2a2b1":["B5a2a2b1a"], "B5b":["B5b1","B5b2","B5b3","B5b4","B5b5"], "B5b1":["B5b1a","B5b1c"], "B5b1a":["B5b1a1","B5b1a2"], "B5b1a2":["B5b1a2a"], "B5b1c":["B5b1c1"], "B5b1c1":["B5b1c1a"], "B5b2":["B5b2a","B5b2b","B5b2c"], "B5b2a":["B5b2a1","B5b2a2"], "B5b2a2":["B5b2a2a"], "B5b2a2a":["B5b2a2a1","B5b2a2a2"], "B5b2c":["B5b2c1"], "B5b3":["B5b3a","B5b3b"]}

            mt_tree_HV = {"HV":["HV0","HV1","HV2","HV20","HV4","HV5","HV6","HV7","HV8","HV9","HV10","HV11","HV14","HV15","HV16","HV17","HV22","HV23","HV24","HV12","HV13","HV18","HV19","HV21", "H"], "HV0":["HV0a","HV0b","HV0c","HV0d","HV0e","HV0f","HV0g"], "HV0a":["HV0a1","V"], "HV0a1":["HV0a1a"], "HV1":["HV1a'b'c","HV1d"], "HV1a'b'c":["HV1a","HV1b","HV1c"], "HV1a":["HV1a1","HV1a2","HV1a3"], "HV1a1":["HV1a1a","HV1a1b"], "HV1a2":["HV1a2a","HV1a2b"], "HV1a3":["HV1a3a"], "HV1b":["HV1b1","HV1b2","HV1b3"], "HV1b1":["HV1b1a","HV1b1b"], "HV1b3":["HV1b3a","HV1b3b"], "HV2":["HV2a"], "HV2a":["HV2a1","HV2a2","HV2a3"], "HV4":["HV4a","HV4b","HV4c"], "HV4a":["HV4a1","HV4a2"], "HV4a1":["HV4a1a"], "HV4a1a":["HV4a1a1","HV4a1a2","HV4a1a3","HV4a1a4"], "HV4a2":["HV4a2a","HV4a2b"], "HV5":["HV5a","HV5b"], "HV6":["HV6a"], "HV9":["HV9a","HV9b","HV9c"], "HV9a":["HV9a1"], "HV9a1":["HV9a1a"], "HV11":["HV11a"], "HV14":["HV14a"], "HV17":["HV17a"], "HV12":["HV12a","HV12b"], "HV12a":["HV12a1"], "HV12b":["HV12b1"], "HV12b1":["HV12b1a"], "HV13":["HV13a","HV13b"]}

            mt_tree_H = {"H":["H1","H2","H3","H4","H5'36","H6","H7","H8","H31","H11","H12","H91","H108","H9","H32","H46","H52","H69","H103","H107","H10","H13","H14","H15","H16","H17","H27","H18","H19","H20","H21","H22","H23","H24","H25","H26","H28","H29","H30","H33","H34","H64","H85","H35","H39","H40","H41","H42","H43","H44","H45","H47","H48","H49","H50","H51","H53","H54","H55","H56","H57","H58","H59","H60","H61","H62","H63","H65","H66","H67","H70","H71","H72","H73","H74","H75","H76","H77","H78","H79","H80","H81","H82","H83","H84","H87","H86","H88","H89","H90","H92","H93","H94","H95","H96","H99","H100","H101","H102","H104","H105","H106"], "H1":["H1a","H1b","H1f","H1g","H1k","H1y","H1z","H1aa","H1ab","H1ac","H1ad","H1cc","H1c","H1e","H1h","H1i","H1an","H1bb","H1j","H1m","H1n","H1o","H1ck","H1p","H1q","H1r","H1s","H1t","H1u","H1v","H1w","H1x","H1ae","H1af","H1ag","H1ah","H1ai","H1aj","H1ak","H1am","H1ao","H1cg","H1ap","H1aq","H1ar","H1as","H1at","H1au","H1av","H1aw","H1ax","H1ay","H1az","H1ba","H1bc","H1bd","H1be","H1bf","H1bg","H1bh","H1ch","H1bj","H1bi","H1bk","H1bm","H1bn","H1bo","H1bp","H1bq","H1br","H1bs","H1bt","H1bu","H1bv","H1bw","H1bx","H1bz","H1ca","H1cd","H1cf","H1ci","H1cj"], "H1a":["H1a1","H1a2","H1a3","H1a4","H1a5","H1a6","H1a7","H1a8","H1a9"], "H1a1":["H1a1a","H1a1b","H1a1c"], "H1a1a":["H1a1a1"], "H1a3":["H1a3a","H1a3b","H1a3c","H1a3d"], "H1a3a":["H1a3a1","H1a3a2","H1a3a3","H1a3a4"], "H1a3b":["H1a3b1"], "H1a3c":["H1a3c1"], "H1a8":["H1a8a"], "H1b":["H1b1","H1b2","H1b3","H1b4","H1b5"], "H1b1":["H1b1a","H1b1b","H1b1c","H1b1d","H1b1h","H1b1e","H1b1f","H1b1g","H1b1h","H1b1i"], "H1b1e":["H1b1e1"], "H1b2":["H1b2a"], "H1b2a":["H1b2a1"], "H1f":["H1f1"],"H1f1":["H1f1a"], "H1g":["H1g1","H1g2"], "H1k":["H1k1","H1k1a"], "H1z":["H1z1"], "H1aa":["H1aa1"], "H1ab":["H1ab1"], "H1c":["H1c1","H1c2","H1c3","H1c4","H1c5","H1c6","H1c7","H1c8","H1c9","H1c10","H1c11","H1c12","H1c13","H1c14","H1c15","H1c16","H1c17","H1c18","H1c19","H1c20","H1c21","H1c22"], "H1c1":["H1c1a","H1c1b","H1c1c","H1c1d"], "H1c1a":["H1c1a1"], "H1c2":["H1c2a"], "H1c3":["H1c3a","H1c3b"], "H1c4":["H1c4a","H1c4b"], "H1c4a":["H1c4a1"], "H1c4b":["H1c4b1"], "H1c5":["H1c5a"], "H1c9":["H1c9a"], "H1e":["H1e1","H1e2","H1e3","H1e4","H1e5","H1e6","H1e7","H1e8"], "H1e1":["H1e1a","H1e1b","H1e1c"], "H1e1a":["H1e1a1","H1e1a2","H1e1a3","H1e1a4","H1e1a5","H1e1a6","H1e1a7","H1e1a8"], "H1e1b":["H1e1b1"], "H1e1b1":["H1e1b1a","H1e1b1b"], "H1e2":["H1e2a","H1e2b","H1e2c","H1e2d"], "H1e4":["H1e4a"], "H1e5":["H1e5a","H1e5b"], "H1e8":["H1e8a"], "H1h":["H1h1","H1h2"], "H1i":["H1i1","H1i2"], "H1i2":["H1i2a"], "H1an":["H1an1","H1an2"], "H1an1":["H1an1a"], "H1j":["H1j1","H1j2","H1j3","H1j4","H1j5","H1j6","H1j7","H1j8","H1j9"], "H1j1":["H1j1a","H1j1b","H1j1c"], "H1j1a":["H1j1a1","H1j1a2"], "H1j2":["H1j2a"], "H1m":["H1m1"], "H1n":["H1n1","H1n2","H1n3","H1n4","H1n5","H1n6"], "H1n1":["H1n1a","H1n1b"], "H1q":["H1q1","H1q2","H1q3"], "H1q1":["H1q1a"], "H1r":["H1r1"], "H1s":["H1s1"], "H1t":["H1t1","H1t2"], "H1t1":["H1t1a"], "H1t1a":["H1t1a1"], "H1u":["H1u1","H1u2"], "H1v":["H1v1"], "H1v1":["H1v1a","H1v1b"], "H1ae":["H1ae1","H1ae2","H1ae3"], "H1ae2":["H1ae2a"], "H1ae3":["H1ae3a"], "H1af":["H1af1","H1af2"], "H1af1":["H1af1a","H1af1b"], "H1ag":["H1ag1"], "H1ag1":["H1ag1a","H1ag1b"], "H1ah":["H1ah1","H1ah2"], "H1ai":["H1ai1"], "H1aj":["H1aj1"], "H1aj1":["H1aj1a"], "H1ak":["H1ak1","H1ak2"], "H1am":["H1am1"], "H1ao":["H1ao1"], "H1ap":["H1ap1"], "H1aq":["H1aq1"], "H1ar":["H1ar1"], "H1as":["H1as1","H1as2"], "H1as1":["H1as1a"], "H1at":["H1at1"], "Htat1":["H1at1a"], "H1au":["H1au1"], "H1au1":["H1au1a","H1au1b"], "H1av":["H1av1"], "H1av1":["H1av1a"], "H1aw":["H1aw1"], "H1ax":["H1ax1"], "H1ba":["H1ba1"], "H1bf":["H1bf1"], "H1bt":["H1bt1"], "H1bv":["H1bv1"], "H2":["H2a","H2b","H2c"], "H2a":["H2a1","H2a2","H2a3","H2a4","H2a5"], "H2a1":["H2a1a","H2a1b","H2a1c","H2a1d","H2a1e","H2a1f","H2a1g","H2a1i","H2a1j","H2a1k","H2a1m","H2a1n"], "H2a1a":["H2a1a1","H2a1a2"], "H2a1b":["H2a1b1","H2a1b2"], "H2a1e":["H2a1e1"], "H2a1e1":["H2a1e1a","H2a1e1b"], "H2a1e1a":["H2a1e1a1"], "H2a1f":["H2a1f1","H2a1f2"], "H2a2":["H2a2a","H2a2b"], "H2a2a":["H2a2a1","H2a2a2"], "H2a2a1":["H2a2a1a","H2a2a1b","H2a2a1c","H2a2a1d","H2a2a1e","H2a2a1f","H2a2a1g","H2a2a1h"], "H2a2b":["H2a2b1","H2a2b2","H2a2b3","H2a2b4","H2a2b5"], "H2a2b1":["H2a2b1a"], "H2a2b1a":["H2a2b1a1"], "H2a2b5":["H2a2b5a"], "H2a3":["H2a3a","H2a3b"], "H2a3a":["H2a3a1"], "H2a5":["H2a5a","H2a5b"], "H2a5a":["H2a5a1"],"H2a5a1":["H2a5a1a","H2a5a1b"], "H2a5b":["H2a5b1","H2a5b2"], "H2c":["H2c1"], "H3":["H3a","H3g","H3i","H3j","H3k", "H3b","H3c","H3d","H3e","H3h","H3m","H3n","H3p","H3q","H3r","H3s","H3t","H3u","H3v","H3w","H3x","H3y","H3z","H3aa","H3ab","H3ac","H3ad","H3ae","H3af","H3ag","H3ah","H3ai","H3aj","H3ak","H3am","H3an","H3ao", "H3ap","H3aq","H3ar","H3as","H3at","H3av"], "H3a":["H3a1"], "H3a1":["H3a1a"], "H3g":["H3g1","H3g2","H3g3","H3g4"], "H3g1":["H3g1a","H3g1b"], "H3i":["H3i1"], "H3k":["H3k1"], "H3k1":["H3k1a"], "H3b":["H3b1","H3b2","H3b3","H3b4","H3b5","H3b6","H3b7"],"H3b1":["H3b1a","H3b1b"], "H3b1b":["H3b1b1"],"H3b1b1":["H3b1b1a"], "H3b4":["H3b4a"], "H3b6":["H3b6a"], "H3c":["H3c1","H3c2","H3c3"], "H3c2":["H3c2a","H3c2b","H3c2c"], "H3c2a":["H3c2a1"], "H3c2b":["H3c2b1"], "H3h":["H3h1","H3h2","H3h3","H3h4","H3h5","H3h6","H3h7"], "H3h2":["H3h2a"], "H3h3":["H3h3a","H3h3b"], "H3q":"H3q1", "H3r":["H3r1"], "H3u":["H3u1"], "H3v":["H3v1","H3v2"], "H3x":["H3x1"], "H3z":["H3z1","H3z2"], "H3ag":["H3ag1"], "H3ao":["H3ao1"], "H3at":["H3at1"], "H4":["H4a","H4b","H4c","H4d"], "H4a":["H4a1","H4a2"], "H4a1":["H4a1a","H4a1c","H4a1d"], "H4a1a":["H4a1a1","H4a1a2","H4a1a3","H4a1a4","H4a1a5"], "H4a1a1":["H4a1a1a"],"H4a1a1a":["H4a1a1a1","H4a1a1a2","H4a1a1a3","H4a1a1a4"], "H4a1a1a1":["H4a1a1a1a"], "H4a1a1a1a":["H4a1a1a1a1"], "H4a1a2":["H4a1a2a"], "H4a1a2a":["H4a1a2a1"], "H4a1a3":["H4a1a3a"], "H4a1a4":["H4a1a4a","H4a1a4b"], "H4a1a4b":["H4a1a4b1","H4a1a4b2"], "H4a1c":["H4a1c1","H4a1c2"], "H4a1c1":["H4a1c1a"], "H4b":["H4b1"], "H4c":["H4c1"], "H5'36":["H5","H36"], "H5":["H5a","H5b","H5c","H5d","H5e","H5f","H5g","H5h","H5j","H5k","H5m","H5n","H5p","H5q","H5r","H5s","H5t","H5u","H5v"], "H5a":["H5a1","H5a2","H5a3","H5a4","H5a5","H5a6","H5a7","H5a8","H5a9"], "H5a1":["H5a1a","H5a1b","H5a1c","H5a1d","H5a1e","H5a1f","H5a1g","H5a1h","H5a1i","H5a1j","H5a1k","H5a1m","H5a1n","H5a1p","H5a1q"], "H5a1c":["H5a1c1","H5a1c2"],"H5a1c1":["H5a1c1a"], "H5a1g":["H5a1g1","H5a1g2"], "H5a1g1":["H5a1g1a"], "H5a2":["H5a2a"], "H5a3":["H5a3a","H5a3b"], "H5a3a":["H5a3a1","H5a3a2","H5a3a3"], "H5a4":["H5a4a"],"H5a4a":["H5a4a1"],"H5a4a1":["H5a4a1a"], "H5a6":["H5a6a"],"H5b":["H5b1","H5b2","H5b3","H5b4","H5b5"], "H5c":["H5c1","H5c2"], "H5c1":["H5c1a"], "H5e":["H5e1"], "H5e1":["H5e1a","H5e1b"], "H5e1a":["H5e1a1"], "H5r":["H5r1","H5r2"], "H5u":["H5u1"], "H6":["H6a","H6b","H6c"], "H6a":["H6a1","H6a2"], "H6a1":["H6a1a","H6a1b"], "H6a1a":["H6a1a1","H6a1a2","H6a1a3","H6a1a4","H6a1a5","H6a1a6","H6a1a7","H6a1a8","H6a1a9","H6a1a10"], "H6a1a1":["H6a1a1a"], "H6a1a2":["H6a1a2a","H6a1a2b"], "H6a1a2b":["H6a1a2b1"], "H6a1a3":["H6a1a3a"], "H6a1a8":["H6a1a8a"], "H6a1b":["H6a1b1","H6a1b2","H6a1b3","H6a1b4"], "H6a1b2":["H6a1b2a","H6a1b2b","H6a1b2c","H6a1b2d","H6a1b2e"], "H6a1b3":["H6a1b3a","H6a1b3b"], "H6a2":["H6a2a"], "H6b":["H6b1","H6b2"], "H6c":["H6c1"], "H7":["H7a","H7b","H7c","H7d","H7e","H7f","H7g","H7h","H7i"], "H7a":["H7a1","H7a2"], "H7a1":["H7a1a","H7a1b","H7a1c","H7a1d"], "H7b":["H7b1","H7b2","H7b3","H7b4","H7b5","H7b6"], "H7b2":["H7b2a"], "H7c":["H7c1","H7c2","H7c3","H7c4","H7c5","H7c6"], "H7d":["H7d1","H7d2","H7d3","H7d4","H7d5"], "H7d2":["H7d2a"], "H7d3":["H7d3a"], "H7h":["H7h1"], "H7i":["H7i1"], "H8":["H8a","H8b","H8c"], "H8a":["H8a1"], "H8b":["H8b1"], "H8c":["H8c1","H8c2"], "H31":["H31a","H31b"], "H11":["H11a","H11b"], "H11a":["H11a1","H11a2","H11a3","H11a4","H11a5","H11a6","H11a7","H11a8"], "H11a2":["H11a2a"], "H11a2a":["H11a2a1","H11a2a2","H11a2a3"], "H11b":["H11b1"], "H12":["H12a"], "H9":["H9a"], "H46":["H46a","H46b"], "H10":["H10a","H10b","H10c","H10d","H10e","H10f","H10g","H10h"], "H10a":["H10a1"], "H10a1":["H10a1a","H10a1b"], "H10a1a":["H10a1a1"], "H10b":["H10b1"], "H10c":["H10c1"], "H10e":["H10e1","H10e2","H10e3"], "H10e1":["H10e1a"], "H10e3":["H10e3a"], "H13":["H13a","H13b","H13c"], "H13a":["H13a1","H13a2"], "H13a1":["H13a1a","H13a1b","H13a1c","H13a1d"], "H13a1a":["H13a1a1","H13a1a2","H13a1a3","H13a1a4","H13a1a5","H13a1a6"], "H13a1a1":["H13a1a1a","H13a1a1b","H13a1a1c","H13a1a1d","H13a1a1e"], "H13a1a1d":["H13a1a1d1"], "H13a1a2":["H13a1a2a","H13a1a2b"], "H13a2":["H13a2a","H13a2b","H13a2c"], "H13a2a":["H13a2a1"], "H13a2b":["H13a2b1","H13a2b2","H13a2b3","H13a2b4","H13a2b5"], "H13a2b2":["H13a2b2a"], "H13a2c":["H13a2c1"], "H13b":["H13b1","H13b2"], "H13b1":["H13b1a","H13b1b"], "H13c":["H13c1","H13c2"], "H13c1":["H13c1a"], "H14":["H14a","H14b"], "H14a":["H14a1","H14a2"], "H14a2":["H14a2a","H14a2b","H14a2c"], "H14b":["H14b1","H14b2","H14b3","H14b4"], "H14b2":["H14b2a"], "H15":["H15a","H15b"], "H15a":["H15a1"], "H15a1":["H15a1a","H15a1b"], "H15a1a":["H15a1a1"], "H15b":["H15b1","H15b2"], "H16":["H16a","H16b","H16c","H16d"],"H16a":["H16a1"], "H17":["H17a","H17b","H17c"], "H17a":["H17a1","H17a2"], "H27":["H27a","H27b","H27c","H27d","H27e","H27f"], "H18":["H18b"], "H20":["H20a","H20b","H20c"], "H20a":["H20a1","H20a2"],"H20a1":["H20a1a"], "H24":["H24a","H24b"], "H24a":["H24a1","H24a2"], "H26":["H26a","H26b","H26c"], "H26a":["H26a1"], "H26a1":["H26a1a","H26a1b"], "H26a1a":["H26a1a1"], "H28":["H28a"], "H28a":["H28a1","H28a2"], "H29":["H29a","H29b"], "H30":["H30a","H30b"], "H30b":["H30b1"], "H33":["H33a","H33b","H33c"], "H35":["H35a"], "H39":["H39a","H39b", "H39c"], "H39a":["H39a1"], "H40":["H40a","H40b"], "H41":["H41a"], "H42":["H42a"], "H42a":["H42a1","H42a2"], "H44":["H44a","H44b"], "H44a":["H44a1"], "H45":["H45a","H45b"], "H47":["H47a"], "H49":["H49a","H49b"], "H49a":["H49a1","H49a2"], "H51":["H51a"], "H55":["H55a","H55b"], "H56":["H56a","H56b","H56c","H56d"], "H56a":["H56a1"], "H58":["H58a"], "H59":["H59a"], "H60":["H60a"], "H61":["H61a"], "H63":["H63a"], "H65":["H65a"], "H66":["H66a"], "H66a":["H66a1"], "H67":["H67a"], "H73":["H73a"], "H73a":["H73a1"], "H76":["H76a"], "H79":["H79a"], "H81":["H81a"], "H95":["H95a"], "H104":["H104a"], "H105":["H105a"]}

            mt_tree_V = {"V":["V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20","V21","V22","V23","V24","V25","V26","V27","V28"], "V1":["V1a","V1b"], "V1a":["V1a1"], "V1a1":["V1a1a","V1a1b"], "V1a1a":["V1a1a1"], "V2":["V2a","V2b","V2c"], "V2a":["V2a1"], "V2a1":["V2a1a"], "V2b":["V2b1","V2b2"], "V3":["V3a","V3b","V3c"], "V3a":["V3a1"], "V7":["V7a","V7b"], "V7a":["V7a1"], "V9":["V9a"], "V9a":["V9a1","V9a2"], "V10":["V10a","V10b"], "V10b":["V10b1","V10b2"], "V15":["V15a"], "V18":["V18a"]}

            mt_tree_F = {"F":["F1","F2","F3","F4"], "F1":["F1a'c'f", "F1b", "F1d","F1e","F1g"], "F1a'c'f":["F1a","F1c","F1f"], "F1a":["F1a1'4","F1a2","F1a3"], "F1a1'4":["F1a1","F1a4"], "F1a1":["F1a1a","F1a1b","F1a1c","F1a1d"], "F1a1a":["F1a1a1"], "F1a1c":["F1a1c1","F1a1c2","F1a1c3"], "F1a1d":["F1a1d1"], "F1a4":["F1a4a","F1a4b"], "F1a4a":["F1a4a1"], "F1a2":["F1a2a"], "F1a3":["F1a3a","F1a3b"], "F1a3a":["F1a3a1","F1a3a2","F1a3a3"], "F1a3a1":["F1a3a1a"], "F1a3a3":["F1a3a3a"], "F1c":["F1c1"], "F1c1":["F1c1a"], "F1c1a":["F1c1a1","F1c1a2"], "F1c1a1":["F1c1a1a","F1c1a1b"], "F1b":["F1b1"], "F1b1":["F1b1a","F1b1b","F1b1c","F1b1d","F1b1e","F1b1f"], "F1b1a":["F1b1a1","F1b1a2"], "F1b1a1":["F1b1a1a"], "F1b1a1a":["F1b1a1a1","F1b1a1a2","F1b1a1a3"], "F1b1a1a1":["F1b1a1a1a"], "F1b1e":["F1b1e1"], "F1d":["F1d1"], "F1e":["F1e1","F1e2","F1e3"], "F1e1":["F1e1a"], "F1g":["F1g1"], "F2":["F2a","F2b","F2g","F2c","F2d","F2e","F2f","F2h","F2i"], "F3":["F3a","F3b"], "F3a":["F3a1"], "F3b":["F3b1"], "F3b1":["F3b1a","F3b1b"], "F3b1a":["F3b1a1","F3b1a2"], "F3b1b":["F3b1b1"], "F4":["F4a","F4b"], "F4a":["F4a1","F4a2"], "F4a1":["F4a1a","F4a1b"], "F4b":["F4b1"]}

            mt_tree_K = {"K":["K1","K2","K3"],
             "K1":["K1a","K1b","K1c","K1d","K1e","K1f"], "K1a":["K1a1","K1a2","K1a3","K1a4","K1a5","K1a6","K1a7","K1a8","K1a9","K1a10","K1a13","K1a14","K1a15","K1a16","K1a26","K1a11","K1a24","K1a30","K1a31","K1a12","K1a17","K1a18","K1a19","K1a23","K1a25","K1a27","K1a28","K1a29"], "K1a1":["K1a1a","K1a1b","K1a1c"], "K1a1a":["K1a1a1","K1a1a2"], "K1a1a2":["K1a1a2a"], "K1a1a2a":["K1a1a2a1"], "K1a1b":["K1a1b1","K1a1b2"], "K1a1b1":["K1a1b1a","K1a1b1b","K1a1b1c","K1a1b1d","K1a1b1e","K1a1b1f","K1a1b1g"], "K1a1b1b":["K1a1b1b1"], "K1a1b2":["K1a1b2a","K1a1b2b"], "K1a1b2a":["K1a1b2a1"], "K1a1b2a1":["K1a1b2a1a"], "K1a2":["K1a2a","K1a2b","K1a2c"], "K1a2a":["K1a2a1","K1a2a2"], "K1a3":["K1a3a"], "K1a3a":["K1a3a1","K1a3a2","K1a3a3","K1a3a4"], "K1a3a1":["K1a3a1a","K1a3a1b"], "K1a4":["K1a4a","K1a4b","K1a4c","K1a4d","K1a4e","K1a4f","K1a4g","K1a4h","K1a4i","K1a4j"], "K1a4a":["K1a4a1"], "K1a4a1":["K1a4a1a","K1a4a1b","K1a4a1c","K1a4a1d","K1a4a1e","K1a4a1f","K1a4a1g","K1a4a1h","K1a4a1i"], "K1a4a1a":["K1a4a1a1","K1a4a1a3","K1a4a1a2"], "K1a4a1a2":["K1a4a1a2a","K1a4a1a2b"], "K1a4a1b":["K1a4a1b1","K1a4a1b2"], "K1a4a1c":["K1a4a1c1"], "K1a4a1f":["K1a4a1f1"], "K1a4b":["K1a4b1"], "K1a4c":["K1a4c1"], "K1a4f":["K1a4f1"], "K1a4h":["K1a4h1"], "K1a4j":["K1a4j1"], "K1a5":["K1a5a","K1a5b"], "K1a8":["K1a8a","K1a8b"], "K1a8a":["K1a8a1"], "K1a10":["K1a10a"], "K1a13":["K1a13a"], "K1a11":["K1a11a","K1a11b"], "K1a11a":["K1a11a1"], "K1a24":["K1a24a"], "K1a30":["K1a30a"], "K1a12":["K1a12a"], "K1a12a":["K1a12a1"], "K1a12a1":["K1a12a1a"], "K1a17":["K1a17a"], "K1a19":["K1a19a"], "K1a29":["K1a29a"], "K1b":["K1b1","K1b2"], "K1b1":["K1b1a","K1b1b","K1b1c"], "K1b1a":["K1b1a1","K1b1a2"], "K1b1a1":["K1b1a1a","K1b1a1b","K1b1a1c","K1b1a1d"], "K1b1a1c":["K1b1a1c1"], "K1b1a1d":["K1b1a1d1"], "K1b1b":["K1b1b1"], "K1b2":["K1b2a","K1b2b"], "K1b2a":["K1b2a1","K1b2a2","K1b2a3"], "K1b2a1":["K1b2a1a"], "K1b2a1a":["K1b2a1a1"], "K1b2a2":["K1b2a2a"], "K1b2b":["K1b2b1"], "K1c":["K1c1","K1c2"], "K1c1":["K1c1a","K1c1b","K1c1c","K1c1d","K1c1e","K1c1f","K1c1g","K1c1h","K1c1i"], "K1c2":["K1c2a"], "K1d":["K1d1"], "K1e":["K1e1"],
             "K2":["K2a","K2b","K2c"], "K2a":["K2a1","K2a2","K2a3","K2a4","K2a5","K2a6","K2a7","K2a8","K2a9","K2a10","K2a11"], "K2a1":["K2a1a"], "K2a2":["K2a2a"], "K2a2a":["K2a2a1"], "K2a3":["K2a3a"], "K2a3a":["K2a3a1"], "K2a5":["K2a5a","K2a5b"], "K2a5a":["K2a5a1"], "K2b":["K2b1","K2b2"], "K2b1":["K2b1a","K2b1b"], "K2b1a":["K2b1a1","K2b1a2","K2b1a3","K2b1a4"], "K2b1a1":["K2b1a1a"]}
            
            try:
                # read the input Eurasian dataset
                haplo_df = pd.read_excel(in_file, header=0)
                
                # repeat the input from the user
                print("\nYou have selected {} on mtDNA".format(haplo_name))
                # search for the individuals belong to the query haplogroup
                query_index = haplo_df[haplo_df['mt_haplogroup']==haplo_name].index.tolist()
                # create a dataframe to contain the selected individuals for later plotting
                query_df = haplo_df.iloc[query_index]
                # check if there is any individual belong to the query haplogroup
                if not query_index:
                    print("There isn't any individual belong to {} in the dataset.".format(haplo_name))
                query = haplo_name
                    
                # by default, the program will only search for three times (if no match in the first closest, it will seach for the second and the third)
                n = 0
                while n<num:
                    # first check if the input haplogroup is in the main trunk (Y-main)
                    if haplo_name in mt_main:
                        haplo_close = get_close(dict_data=mt_tree_main, haplogroup=haplo_name)
                        print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        if haplo_close == "mt-MRCA":
                            n = num
                            print("\nThe search reached the end: Eve mt-MRCA! Please try another haplogroup.")
                        else:
                            # create a list to contain the indice for the matched individuals
                            individual = haplo_df[haplo_df["mt_haplogroup"]==haplo_close].index.tolist()
                        
                            # if there are individuals belong to the closest haplogroup
                            if individual:
                                n = num
                                match_df = haplo_df.iloc[individual]
                                individual_count = len(match_df)
                                print("Found {} individual(s) in the haplogroup {}".format(individual_count,haplo_close))
                                
                                # create the output file name
                                out_file = "mt_" + haplo_close + ".pdf"
                                
                                # plot the results on the map
                                map_plot(out_file, query_df, match_df, query, haplo_close, colors)
                                print("\nThank you for using HaploMap! The graph has been printed to your working directory.")
                                
                            else: # there isn't any individual in the closest group
                                n += 1
                                if n<num:
                                    print("There isn't any individual in the closest haplogroup, searching for next closest!")    
                                    haplo_name = haplo_close
                                else:
                                    print('\nThank you for using HaploMap! We only search for the top {} closest haplogroups. No matched individuals!'.format(num))
                    
                    # the haplogroup is not in the main tree trunk:
                    else:
                        group = re.findall('[A-Z]{1,2}', haplo_name)[0]
                        # go to the corresponding sub-tree for searching (there are totally 27 sub-trees)
                        if group=='A':
                            haplo_close = get_close(dict_data=mt_tree_A, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='B':
                            haplo_close = get_close(dict_data=mt_tree_B, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='C':
                            haplo_close = get_close(dict_data=mt_tree_C, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='D':
                            haplo_close = get_close(dict_data=mt_tree_D, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='E':
                            haplo_close = get_close(dict_data=mt_tree_E, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='F':
                            haplo_close = get_close(dict_data=mt_tree_F, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='G':
                            haplo_close = get_close(dict_data=mt_tree_G, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='HV':
                            haplo_close = get_close(dict_data=mt_tree_HV, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='H':
                            haplo_close = get_close(dict_data=mt_tree_H, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='I':
                            haplo_close = get_close(dict_data=mt_tree_I, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='J':
                            haplo_close = get_close(dict_data=mt_tree_J, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='K':
                            haplo_close = get_close(dict_data=mt_tree_K, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='L':
                            haplo_close = get_close(dict_data=mt_tree_L, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='M':
                            haplo_close = get_close(dict_data=mt_tree_M, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='N':
                            haplo_close = get_close(dict_data=mt_tree_N, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='O':
                            haplo_close = get_close(dict_data=mt_tree_O, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='P':
                            haplo_close = get_close(dict_data=mt_tree_P, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='Q':
                            haplo_close = get_close(dict_data=mt_tree_Q, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='R':
                            haplo_close = get_close(dict_data=mt_tree_R, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='S':
                            haplo_close = get_close(dict_data=mt_tree_S, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='T':
                            haplo_close = get_close(dict_data=mt_tree_T, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='U':
                            haplo_close = get_close(dict_data=mt_tree_U, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='V':
                            haplo_close = get_close(dict_data=mt_tree_V, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='W':
                            haplo_close = get_close(dict_data=mt_tree_W, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='X':
                            haplo_close = get_close(dict_data=mt_tree_X, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='Y':
                            haplo_close = get_close(dict_data=mt_tree_Y, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        elif group=='Z':
                            haplo_close = get_close(dict_data=mt_tree_Z, haplogroup=haplo_name)
                            print("The closest haplogroup should be {}\nsearching for the individuals...".format(haplo_close))
                        
                        # create a list to contain the indice for the matched individuals
                        individual = haplo_df[haplo_df["mt_haplogroup"]==haplo_close].index.tolist()
                    
                        # if there are individuals belong to the closest haplogroup
                        if individual:
                            n = num
                            match_df = haplo_df.iloc[individual]
                            individual_count = len(match_df)
                            print("Found {} individual(s) in the haplogroup {}".format(individual_count,haplo_close))
                            
                            # create the output file name
                            out_file = "mt_" + haplo_close + ".pdf"
                            
                            # plot the results on the map
                            map_plot(out_file, query_df, match_df, query, haplo_close, colors)
                            print("\nThank you for using HaploMap! The graph has been printed to your working directory.")
                            
                        else: # there isn't any individual in the closest group
                            n += 1
                            if n<num:
                                print("There isn't any individual in the closest haplogroup, searching for next closest!")    
                                haplo_name = haplo_close
                            else:
                                print('\nThank you for using HaploMap! We only search for the top {} closest haplogroups. No matched individuals!'.format(num))
            
            except TypeError:
                print('TypeError: The input is not a defined haplogroup on mtDNA. Please check again!')
            
            except IndexError:
                print("The haplogroup input was in wrong format! Please check again!")
            
            
            except FileNotFoundError as not_found:
                print("The file {} was not found!".format(not_found.filename))

        else:
            print("Please enter a correct chromosome! Y or mt!")


        """
2. The second function: input a mutation name, search for its subgroup and plot on the map (only available for Y-DNA now).

Steps:
    1) Read the SNP_index.xlsx as the reference database. Check if the user's input is within the index.
    2) If the mutation is found in the index, extract its subgroup, build 37 number, build 28 number and mutation information from the database.
    3) Read the test_Eurasian.xlsx dataset, seach for the individuals belonging to the subgroup found in step 2.
    4) Plot the individuals on the map.

Example usage: 
    -for Y-DNA: python HaploMap.py --mode 2 --input test_Eurasian.xlsx 
        Please enter the mutation name: V1023

        """
    
    elif mode == '2':
        mut_name = input("Please enter the mutation name (Y-DNA):")
        
        # read the SNP_index file as a reference database
        SNP_dir = sys.argv[0][:-11] # by default, it will locate at the same directory as the HaploMap.py
        SNP_file = SNP_dir + 'SNP_index.xlsx'
        SNP_df = pd.read_excel(SNP_file, header=0, index_col=0)
        
        # define a color list for plotting
        intervals = ['1001-2000 CE', '1-1000 CE', '1000-1 BCE', '2000-1001 BCE', '3000-2001 BCE', '4000-3001 BCE', '5000-4001 BCE', '6000-5001 BCE', '7000-6001 BCE', '8000-7001 BCE', '9000-8001 BCE', '10000-9001 BCE', '11000-10001 BCE']
        color_list = ['lightcoral','brown','red','darkorange','gold','yellowgreen','limegreen','blue','violet','fuchsia','darkorchid','yellow','cyan']
        colors = {}
        for i in range(len(intervals)):
            key = intervals[i]
            colors[key] = color_list[i]
        
        try:
            
            # to check if the mutation is included in the SNP_index
            if mut_name in SNP_df.index.values:
                # extract the haplogroup from the SNP_index
                haplogroup = SNP_df.loc[mut_name,'Subgroup Name']
                build37 = str(SNP_df.loc[mut_name, 'Build 37 Number']) # the original class would be int
                build38 = str(SNP_df.loc[mut_name, 'Build 38 Number']) # the original class would be int, and some missing value would be a float
                mut_info = str(SNP_df.loc[mut_name, 'Mutation Info']) # for some missing value, it would be a float
                
                # to adjust the output text for some missing values
                if build38 == 'nan':
                    build38 = "None"
                if mut_info == 'nan':
                    mut_info = "Not specified"
                # print the information to the screen to announce the user
                print("\nQuery mutation name: {}\nHaplogroup Name: {}\nGRCh37 (Build 37 number): {}\nGRCh38 (Build 38 number): {}\nMutation information: {}\n\nNow searching for the individuals in the {} dataset...".format(mut_name, haplogroup, build37, build38, mut_info, in_file))
                
                # read the dataset
                haplo_df = pd.read_excel(in_file, header=0)
                # create a match_df dataframe to extract the individuals belong to the haplogroup
                match_df = haplo_df.loc[haplo_df['Y_haplogroup'] == haplogroup]
                individual_count = len(match_df)
                print("\nFound {} individual(s) in the {} dataset.".format(individual_count, in_file))
                
                # specify the output report's file name
                out_report = mut_name + '.report.txt'
                
                # if at least one individual were found in the haplogroup
                if individual_count>0:
                    # specify the output plot's file name
                    out_plot = mut_name + '.pdf'
                    
                    # plot the individuals in the defined haplogroup on the map
                    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
                    fig,ax=plt.subplots(dpi=600)
                    ax.set_aspect('equal')
                    ax.xaxis.label.set_visible(False)
                    ax.yaxis.label.set_visible(False)
                    world.plot(ax=ax, color = 'lightgrey', edgecolor = 'darkgrey', linewidth = 0.5)
                    grouped = match_df.groupby('Age_interval')
                    for key, group in grouped:
                        label_name = haplogroup + " (" + key + ")"
                        group.plot(ax=ax, kind='scatter', x='Long.', y='Lat.', s=3.5, label=label_name, color=colors[key])
                    plt.legend(bbox_to_anchor=(1.0,1.0), loc=2, fontsize=8)
                    plt.savefig(out_plot, format="pdf", bbox_inches='tight')
                    
                    # print the results to the report
                    with open(out_report, 'w') as out_report:
                        print("\nQuery mutation name: {}\nHaplogroup Name: {}\nGRCh37 (Build 37 number): {}\nGRCh38 (Build 38 number): {}\nMutation information: {}\n".format(mut_name, haplogroup, build37, build38, mut_info), file = out_report)
                        print("Found {} individual(s) in the {} dataset.".format(individual_count, in_file), file = out_report)
                    print("\nThank you for using HaploMap! Your results are available in the working directory.")
                
                else:
                    # only print the report without a map
                    with open(out_report, 'w') as out_report:
                        print("\nQuery mutation name: {}\nHaplogroup Name: {}\nGRCh37 (Build 37 number): {}\nGRCh38 (Build 38 number): {}\nMutation information: {}\n".format(mut_name, haplogroup, build37, build38, mut_info), file = out_report)
                        print("Found {} individual(s) in the {} dataset.".format(individual_count, in_file), file = out_report)
                    print("\nThank you for using HaploMap! Your report is available in the working directory.")
                
                
            else:
                print("Please check the input of the mutation again! It is not included in the SNP_index.\nIt may due to the wrong mutation name or it was not a mutation mark to define a subgroup.")
            


        except FileNotFoundError as not_found:
            print("The file {} was not found!".format(not_found.filename))



        
        
        """
3. The third function: input a country, calculate the main haplogroup frequency.

Steps:
    1) Define a function haplo_freq() to select the individuals from the selected country in the dataframe, extract their main haplogroups, and store the information into a group_list for easier calculation. The group_list was turned to a set group_uniq to remove repeated groups. The calculation of the frequency was then performed based on the group_set and group_list.
    2) Seperate the output information based on the chromosome chosed by the user.

Example usage: 
    -for Y-DNA: python HaploMap.py --mode 3 --input test_Eurasian.xlsx 
        Please select the chrmosome (Y/mt): Y
        Please select a country to discover: China
    
    -for mt-DNA: python HaploMap.py --mode 3 --input test_Eurasian.xlsx 
        Please select the chrmosome (Y/mt): mt
        Please select a country to discover: China

        """       
        
        
        
    elif mode == '3':
        chr_name = input("Please select the chromosome (Y/mt):")
        country_name = input("Please select a country to discover:")
        
        
        def haplo_freq(df, chr_name, country_name):
            individual = df[df["Country"]==country_name].index.tolist()
            chr_col = chr_name + "_haplogroup"
            match_df = df.loc[individual,['Country', chr_col]]
            match_df = match_df.dropna(axis=0, subset=[chr_col])
            
            # create a new column to contain the main haplogroup
            match_df["Main"] = np.nan
            
            # extract the indice 
            match_index = match_df.index.tolist()
            
            # enter the main haplogroup for each individual
            for i in match_index:
                chr_group = match_df.loc[i, chr_col]
                match_df.loc[i, "Main"] = re.findall("^[A-Z]+", chr_group)[0]
            
            # extract the groups
            group_list = match_df["Main"].tolist()
            group_uniq = set(group_list)
            total_pop = len(group_list)
            
            # calculate and print the results to the output file
            out_file = country_name + "." + chr_name + ".txt"
            with open(out_file, 'w') as out_file:
                for group in group_uniq:
                    count = group_list.count(group)
                    freq = count / total_pop * 100
                    freq = round(freq, 2) # keep two decimals
                    print("{}: {} ({}%)".format(group, count, freq), file=out_file)
                print("total individuals (with {}-DNA information): {}".format(chr_name, total_pop), file=out_file)
            
            
        try:
            haplo_df = pd.read_excel(in_file, header=0)
            
            # check if the input for country is correct.
            if (haplo_df["Country"].eq(country_name)).any():
                
                # for Y-DNA
                if chr_name.upper() == "Y":
                    chr_name = "Y"
                    print("\nCalculating haplogroup frequency on Y-chromosome in {}...".format(country_name))
                    haplo_freq(haplo_df, chr_name, country_name)
                    print("\nThank you for using HaploMap! Your report is available in the working directory.")
                
                # for mt_DNA
                elif chr_name.lower() == "mt":
                    chr_name = "mt"
                    print("\nCalculating haplogroup frequency on mt-DNA in {}...".format(country_name))
                    haplo_freq(haplo_df, chr_name, country_name)
                    print("\nThank you for using HaploMap! Your report is available in the working directory.")
                
                else:
                    raise TypeError
                    
            else:
                print("The country is not included in the {} dataset. Please check the coutry name!".format(in_file))
        
        
        except TypeError:
            print("Please enter a right chromosome (Y/mt)! Or enter --help to check the usage!")
        
        except FileNotFoundError as not_found:
            print("The file {} was not found!".format(not_found.filename))
        
    
    else:
        print("Please select the correct mode! Enter --help to check the usage!")
        
else:
    print("Argument missing! Please enter '--help' to check the usage!")
