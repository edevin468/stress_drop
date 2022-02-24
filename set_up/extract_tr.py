#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: emmadevin

***this script must be run inside gmprocess in order to use the StreamWorkspace module***

inputs: *.h5 files for each event from gmprocess that contain processed waveforms and tells whether or not waveforms passed QA

    saves all waveforms that passed QA as text files (*.out) for each channel into event folders

outputs: *.out text files in the processed directory

"""

import matplotlib.pyplot as plt
plt.style.use("ggplot")
from obspy import read

import os
import os.path as path
import glob
import numpy as np
from spec_func import bin_spec
from spec_func import bin_max_err
import time
from gmprocess.io.asdf.stream_workspace import StreamWorkspace

#working directory here
working = '/Users/emmadevin/Work/USGS_2021/Data/gmprocess/qa_processing_test8/data'
outpath = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/processed'

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
            
            
            tr_E = 0
            tr_N = 0
            tr_Z = 0
            for tr in st:
                channel = tr.stats.channel[-1]
                if channel == 'E': 
                    tr_E = tr
                    data_E = tr_E.data
                    delta_E = 1/tr_E.stats.sampling_rate
                    
                    outfile = open(outpath + '/' +  ev + '/' + stn + '_'  + str(delta_E) + '_E.out', 'w')
                    data_E = data_E.T
                    np.savetxt(outfile, data_E, fmt=['%E'], delimiter='\t')
                    outfile.close()
                    
                
                    
                    
                elif channel == 'N': 
                    tr_N = tr
                    data_N = tr_N.data
                    delta_N = 1/tr_N.stats.sampling_rate
                    
                    outfile = open(outpath + '/' +  ev + '/' + stn + '_'  + str(delta_N) + '_N.out', 'w')
                    data_N = data_N.T
                    np.savetxt(outfile, data_N, fmt=['%E'], delimiter='\t')
                    outfile.close()                
                
                    
                elif channel == 'Z': 
                    tr_Z = tr
                    data_Z = tr_Z.data
                    delta_Z = 1/tr_Z.stats.sampling_rate
                    
                    outfile = open(outpath + '/' +  ev + '/' + stn + '_' + str(delta_Z) + '_Z.out', 'w')
                    data_Z = data_Z.T
                    np.savetxt(outfile, data_Z, fmt=['%E'], delimiter='\t')
                    outfile.close()                
                

                
              