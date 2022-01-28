#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 11:59:01 2022

@author: emmadevin
"""

import os
import pkg_resources
import os.path as path
import glob
from gmprocess.io.asdf.stream_workspace import StreamWorkspace

working = '/Users/emmadevin/Work/USGS_2021/Data/gmprocess/qa_processing_test2/data'
outpath = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/RC_beta'
events = glob.glob(working + '/*')

for event in events:
    name = path.basename(event)
    ev = name.split('/')[-1]
    print(ev)

    data_path =  working + '/' + ev + '/workspace.h5'
    label = ev + '_default'
    
    if not path.exists(outpath + '/' + ev ):
        os.makedirs(outpath + '/' + ev )
        

    workspace = StreamWorkspace.open(data_path)
    
    ds = workspace.dataset
    stn_list = ds.waveforms.list()
    
    for stn in stn_list:
        
        sc = workspace.getStreams(ev, stations=[stn], labels=['default'])
      
        
        sta_st = sc[0]
        check = sta_st.passed
        
        if check == True: 
            st = ds.waveforms[stn][label]
    
            for tr in st:
                filename = stn+'_'+tr.stats.channel+'_'+ev+'.mseed'
                
                tr.write(outpath + '/' + ev + '/' + filename)
        
