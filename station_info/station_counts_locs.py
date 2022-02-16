#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 10:13:32 2022

@author: emmadevin
"""

import obspy as op
from obspy import read, Stream
import os
import os.path as path
import glob
import pandas as pd


working_dir = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/event_files'

files = glob.glob(working_dir + '/*')
 
dfs = []                 
for file in files:
    
    df = pd.read_csv(file, sep = '|')
    dfs.append(df)
    
event_file = pd.concat(dfs)
    
event_file = event_file[['#Network', 'Station', 'Latitude', 'Longitude']]
event_file = event_file.drop_duplicates(subset=['#Network', 'Station'])
event_file = event_file.reset_index()
event_file = event_file[['#Network', 'Station', 'Latitude', 'Longitude']]
event_file = event_file.drop([307])

file_dir = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/processed'

files  = glob.glob(file_dir + '/*/*')

stn_ids=[]
for file in files: 

    filename = path.basename(file)
      
    stn_id = filename.split('_')[0]
    stn_ids.append(stn_id)

stns =  list(set(stn_ids))

counts = []
sts = []
nets = []
for stn in stns:
    count = stn_ids.count(stn)
    counts.append(count)
    
    st = stn.split('.')[1]
    sts.append(st)  
    net = stn.split('.')[0]
    nets.append(net)
    
df = pd.DataFrame()
df['#Network'] = nets
df['Station']  = sts
df['Count'] = counts 


dff = pd.merge(df, event_file, how = 'inner', on = ['#Network','Station'])

dff.to_csv('/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/station_data/station_counts.csv')






    
