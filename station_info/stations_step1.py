#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 09:19:54 2021

Create *.csv file in each event folder with list of stations 

@author: emmadevin
"""

import obspy as op
from obspy import read, Stream
import os
import os.path as path
import glob
import pandas as pd

# working directory
working_dir = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered'

# event directories and outpath
event_dirs = glob.glob(working_dir + '/RC_beta/*')

# create list of event directory names
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
    
# loop through *.ms files and apply instrument corrections and save the resulting files   
stn_ids = []      
for event in events: 
    
    # event directory
    event_dir = working_dir + '/RC_beta/' + event
    
    # create list of all files in event directory
    file_list = glob.glob(event_dir + '/*.mseed')
    



    for file in file_list: 
        # read in file and determine station 
        st = read(file)
        filename = path.basename(file)
      
        stn_id = filename.split('_')[0]
        stn_ids.append(stn_id)
        
stns =  list(set(stn_ids))
counts = []
for stn in stns:
    count = stn_ids.count(stn)
    counts.append(count)
            
    

df = pd.DataFrame()
df['station_id'] = stns
df['count'] = counts
            
df.to_csv(working_dir + '/station_counts.csv')