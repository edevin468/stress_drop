#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 10:35:04 2022

@author: emmadevin
"""

import obspy as op
from obspy import read
import os
import os.path as path
import glob
import pandas as pd



# working directory
working_dir = '/Users/emmadevin/Work/USGS_2021/Data/gmprocess/event_downloads/data'

# event directories
event_dirs = glob.glob(working_dir + '/*')

outpath = '/Users/emmadevin/Work/USGS_2021/Data/gmprocess/qa_processing_test2/data'


 
for event in event_dirs:
    name = path.basename(event)
    ev = name.split('/')[-1]
    print(ev)
    
    if not path.exists(outpath + '/' + ev + '/raw'):
        os.makedirs(outpath + '/'  + ev + '/raw')
        
    source = working_dir + '/' + ev + '/*.json'
    dest = outpath + '/'+ ev 
    
    copy = 'cp -R ' + source + ' ' + dest
    print(copy)
    
    os.popen('cp -R ' + source + ' ' + dest)
    
    
