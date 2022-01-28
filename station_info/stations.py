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
working_dir = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered'

# event directories and outpath
event_dirs = glob.glob(working_dir + '/RC_beta/*')


# create list of event directory names
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
 
df = pd.read_csv(working_dir + '/RC_beta/ci38443183/38443183.iris.txt', sep = '|')

# loop through events and check which station appears in which response file      
for event in events: 
    
    # event directory
    event_dir = working_dir + '/RC_beta/' + event
    
    # create list of all files in event directory
    file_list = glob.glob(event_dir + '/*.txt')

    for file in file_list:
        inv = pd.read_csv(file, sep = '|')
        
        df = pd.concat([df,inv])
        
network = df['#Network']
station = df['Station']
lat = df['Latitude']
lon = df['Longitude']
        
station_locs = pd.DataFrame()
station_locs['network'] = network
station_locs['station'] = station
station_locs['latitude'] = lat
station_locs['longitude'] =lon

sl = station_locs.drop_duplicates()
sl = sl.drop([0])

# sl.to_csv('/Users/emmadevin/Work/USGS 2021/Data/Prelim+/Station_info/station_locs.csv')