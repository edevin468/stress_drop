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
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim+'

# event directories
event_dirs = glob.glob(working_dir + '/RC_beta/*')
outpath = '/Users/emmadevin/Work/USGS 2021/Data/gmprocess/qa_processing_test2/data'


 
for event in event_dirs:
    name = path.basename(event)
    ev = name.split('/')[-1]
    print(ev)

    os.makedirs(outpath1 + '/ci'  + ev)
    
    
