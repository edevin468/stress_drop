#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 10:24:22 2021

@author: emmadevin
"""

import obspy as op
from obspy import read
import os
import os.path as path
import glob
import pandas as pd



working_dir =  '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered'

events = glob.glob(working_dir + '/Andrews_inversion/Events/*')


id_list = []
station_count_list = []
lat_list = []
lon_list = []
depth_list = []
year_list = []
month_list = []
day_list = []
hour_list = []
min_list = []
sec_list = []
M_list = []

for event in events: 
    
    # get event id
    filename = (event.split('/')[-1])
    event_id = filename.split('.')[0]
    id_list.append(event_id)
    
    # get station count
    file_list = glob.glob(working_dir +'/corrected/' + event_id + '/*')
    station_list = []
    for file in file_list:
        name = file.split('/')[-1]
        ntwk = file.split('.')[0]
        stn = file.split('.')[1]
        stn_id = ntwk + stn
        
        if stn_id not in station_list: station_list.append(stn_id)
        
    station_count = len(station_list)
    station_count_list.append(station_count)
    
    # get linfo from phase file
    phase_file = working_dir + '/RC_phase_beta/' + event_id + '.phase'
    phase = pd.read_csv(phase_file, sep = '\s+', index_col=0, nrows = 0).columns.tolist()
    
    lat = float(phase[4])
    lat_list.append(lat)
    
    lon = float(phase[5])
    lon_list.append(lon)
    
    depth = float(phase[6])
    depth_list.append(depth)
    
    M = float(phase[7])
    M_list.append(M)
   
    date_string = phase[3]
    
    date = date_string.split(',')[0]
    time = date_string.split(',')[1]
    
    year = date.split('/')[0]
    year_list.append(year)
    
    month = date.split('/')[1]
    month_list.append(month)
    
    day = date.split('/')[2]
    day_list.append(day)
    
    hour = time.split(':')[0]
    hour_list.append(hour)
    
    minute = time.split(':')[1]
    min_list.append(minute)
    
    sec = time.split(':')[2]
    sec_list.append(sec)
    


df = pd.DataFrame()

df['event'] = id_list
df['latitude'] = lat_list
df['longitude'] = lon_list
df['depth'] = depth_list
df['year'] = year_list
df['month'] = month_list
df['day'] = day_list
df['hour'] = hour_list
df['min'] = min_list
df['sec'] = sec_list
df['catalog M'] = M_list
df['station_count'] = station_count_list

df.to_csv(working_dir + '/stress_drops/event_info.csv')



    
    