#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 09:26:40 2021

add columns to *.csv files created in stations_step1.py for whether or not the station data exists in each of the *.txt files

@author: emmadevin
"""

import obspy as op
from obspy import read, Stream
import os
import os.path as path
import glob
import pandas as pd

# START FUNCTION DEFINITIONS #
# ==================== #

def search(file_name, string):
    """ Check if any line in the file contains given string """
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            if string in line:
                return True
    return False

# ==================== #
# END FUNCTION DEFINITIONS #

# working directory
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim_qa_filtered'

# event directories and outpath
event_dirs = glob.glob(working_dir + '/RC_beta/*')


# create list of event directory names
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
    
# loop through events and check which station appears in which response file      
for event in events: 
    
    # event directory
    event_dir = working_dir + '/RC_beta/' + event
    
    # create list of all files in event directory
    file_list = glob.glob(event_dir + '/*.txt')
    df = pd.read_csv(event_dir + '/stations.csv')
    stns = df['station_id'].tolist()
    

    for file in file_list:
        check = []
        filename = path.basename(file)
        typ = filename.split('.')[1]
        
        for stn in stns:
            if search(file,str(stn)):
                check.append(True)
            else:
                check.append(False)
                
        df[str(typ)] = check
        
        

    type_list = []
    for i in range(len(df)):
        scedc = df['scedc'].tolist()
        ncedc = df['ncedc'].tolist()
        iris = df['iris'].tolist()
        
        if scedc[i]==True: type_list.append('scedc')
        elif ncedc[i]==True: type_list.append('ncedc')
        elif iris[i]==True: type_list.append('iris')
        else: type_list.append('none')
        
    df['type'] = type_list
        
    


 
    df.to_csv('/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered/RC_beta/'+event+'/station_inv.csv')
    
    
    
    