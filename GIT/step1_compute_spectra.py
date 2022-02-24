#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: emmadevin

Originally written by Alexis Klimasewski.
Modified for use in Ridgecrest stress drop project by Emma Devin:
- uses both N and E channels 
- takes *.out files from gmprocess
- channel combination changed to

***this script must be run outside of spyder due to compatability issue with conda and mtspec***

input: *.out files from gmprocess

    calls bin_spec function to return evenly spaced log bins and binned data,
    takes the average of the N and E components

outputs: writes bins and binned spectra into the record_spectra directory

"""

import matplotlib.pyplot as plt
plt.style.use("ggplot")
from obspy import read
from mtspec import mtspec
import os
import os.path as path
import glob
import numpy as np
from spec_func import bin_spec
from spec_func import bin_max_err
import time

working_dir = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered'

event_dirs = glob.glob(working_dir + '/processed/*')
outpath = working_dir + '/record_spectra'

events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
for i in range(len(events)):
    if not path.exists(outpath + '/' + events[i]):
        os.makedirs(outpath + '/'  + events[i])

for event in events:
    t1 = time.time()
    print('binning and fft of event: ' + event)
    
    recordpaths = glob.glob(working_dir + '/processed/' + event + '/*_N.out')#full path for only specified channel
    stns = [(x.split('/')[-1]).split('_')[0] for x in recordpaths]
    
    for stn in stns:
        recordpath_E = glob.glob(working_dir  + '/processed/' + event +'/' + stn + '*_E.out')
        recordpath_N = glob.glob(working_dir  + '/processed/' + event +'/' + stn + '*_N.out')
        recordpath_Z = glob.glob(working_dir  + '/processed/' + event +'/' + stn + '*_Z.out')
        
        if(len(recordpath_E) == 1 and len(recordpath_N) == 1):
            
            base_N = path.basename(recordpath_N[0])
            base_E = path.basename(recordpath_E[0])
            
            stn_id_N = base_N.split('_')[0]
            network = stn_id_N.split('.')[0]
            station = stn_id_N.split('.')[1]
            
             # spectra and std dev for E channel
            data_E = np.genfromtxt(recordpath_E[0], dtype = float, comments = '#', delimiter = None)
            delta_E = float(base_E.split('_')[1])
            
            
            spec_amp_E, freq, jack_E, fstat_E, dof_E =  mtspec(data_E, delta = delta_E, time_bandwidth = 4, number_of_tapers=7, quadratic = True, statistics = True)
            sigmaE = (jack_E[:,1] - jack_E[:,0])/3.29
            
            spec_array_E = np.array(spec_amp_E)
            freq_array_E = np.array(freq)
            
            # spectra and std dev for N channel
            data_N = np.genfromtxt(recordpath_N[0], dtype = float, comments = '#', delimiter = None)
            delta_N = float(base_N.split('_')[1])
            
            spec_amp_N, freq, jack_N, fstat_N, dof_N =  mtspec(data_N, delta = delta_N, time_bandwidth = 4, number_of_tapers=7, quadratic = True, statistics = True)
            sigmaN = (jack_N[:,1] - jack_N[:,0])/3.29
            
            spec_array_N = np.array(spec_amp_N)
            freq_array_N = np.array(freq)
            
            # if evenly sampled
            if(len(spec_array_E)==len(spec_array_N)):
    
                # # here we bin into evenly spaced bins with frequency
                
                
                # #spectra is power spectra so add the two components
                # data_NE_2 = spec_array_E + spec_array_N
                # #now data is NE power spectra
                # #take the square root for normal velocity spectra
                # data_NE = np.sqrt(data_NE_2)
                
                # log average of components
                data_NE_2 = np.exp((np.log(spec_array_E)+np.log(spec_array_N))/2)
                data_NE = np.sqrt(data_NE_2)
                
                

                sigma = np.sqrt((spec_array_N/data_NE**2.)*sigmaN + ((spec_array_N/data_NE**2.)*sigmaN))
                

                bins, binned_data = bin_spec(data_NE[6:-1], freq[6:-1], num_bins = 75)
                bins_sig, binned_sig = bin_max_err(sigma[6:-1], freq[6:-1], num_bins = 75)

              
                if (np.isnan(binned_data).any() == False):
             
                    outfile = open(outpath + '/' +  event + '/' + network + '_' + station + '_' + 'NE' + '__' + event + '.out', 'w')
                    data = np.array([bins, binned_data, binned_sig])
                    data = data.T
                    outfile.write('#bins \t \t vel_spec_NE_m \t binned_sig \n')
                    np.savetxt(outfile, data, fmt=['%E', '%E', '%E'], delimiter='\t')
                    outfile.close()
                else: print('No data')
                
    t2 = time.time()
    print('time for event: (s)', (t2-t1))
        
