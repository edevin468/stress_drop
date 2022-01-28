#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 10:47:49 2021

removes records from stations with less than 3 records

@author: emmadevin
"""

import obspy as op
from obspy import read, Stream
import os
import os.path as path
import glob
import pandas as pd
import numpy as np


# working directory
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered'

event_dirs = glob.glob(working_dir + '/RC_Beta/*')

ev = glob.glob(working_dir + '/RC_Beta/*/*')

counts = pd.read_csv(working_dir + '/Station_info/station_counts.csv') 

remove_list = []
for event in event_dirs:
    
    
    r = pd.read_csv(event + '/station_inv.csv')
    types = r['type']
    stn_ids = r['station_id']
    
    for i in range(len(types)):
        if types[i] == 'none':
            
            ntwk = stn_ids[i].split('|')[0]
            stn = stn_ids[i].split('|')[1]
            
            remove_list.append(ntwk+'.'+stn)
            
            
counts_list = counts['count']

remove = counts_list<3

stn = counts['station']
ntwk = counts['network']

stn_ids = np.array(ntwk+'.'+stn)

remove_stns = list(stn_ids[remove])

remove_list = remove_list + remove_stns


for event in event_dirs: 
    for stn_id in remove_list:
        files = glob.glob(event + '/' + stn_id + '.*')
        for file in files: 
            if os.path.isfile(file):
                print('True')
                os.remove(file)


