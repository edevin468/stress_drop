#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 13:15:54 2021

@author: emmadevin
"""

import glob
import os.path as path
import numpy as np
import obspy
from obspy import read    
import pandas as pd



working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim+/record_spectra'

file_list = glob.glob(working_dir +'/*/*')
station_list = []
record_list = []

for file in file_list:
    record = (file.split('/')[-1])
    name = path.basename(record)

    stn = name.split('_')[1]
    
    if stn not in station_list:  station_list.append(stn)
    
    record_list.append(stn)


counts = [0]*len(station_list)

for record in record_list:
    
    index = station_list.index(record)
    counts[index] = counts[index] + 1
    
df = pd.DataFrame()
df['station'] = station_list
df['count'] = counts

df.to_csv('/Users/emmadevin/Work/USGS 2021/Data/Prelim+/station_counts.csv')