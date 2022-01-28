#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:06:20 2017

@author: Alexis Klimasewski
Edited by Emma Devin to read .ms files and work with the data we have for Ridgecrest


inputs: reads in paths to the instrument corrected  and cut HHN and HHE channel sac files

method: reads in trace data from sac files
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

#working directory here
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim_qa_filtered'

#path to corrected seismograms
event_dirs = glob.glob(working_dir + '/RC_beta/*')
outpath = working_dir + '/record_spectra'



##make event directories within corrected local data
##make a directory for each event
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
for i in range(len(events)):
    if not path.exists(outpath + '/' + events[i]):
        os.makedirs(outpath + '/'  + events[i])

for event in events:
    t1 = time.time()
    print('binning and fft of event: ' + event)
    
    recordpaths = glob.glob(working_dir + '/RC_beta/' + event + '/*.HHN.*.mseed')#full path for only specified channel
    stns = [(x.split('/')[-1]).split('.')[0] for x in recordpaths]
    
    for stn in stns:
        recordpath_E = glob.glob(working_dir + '/corrected/' + event +'/' + stn + '.*.HHE*.mseed')
        recordpath_N = glob.glob(working_dir + '/corrected/' + event +'/' + stn + '.*.HHN*.mseed')
        
        if(len(recordpath_E) == 1 and len(recordpath_N) == 1):
            #North component
            base_N = path.basename(recordpath_E[0])
            base_E = path.basename(recordpath_N[0])
    
            network = base_N.split('.')[1]
            station = base_N.split('.')[0]
            full_channel_N = base_N.split('.')[2]
            full_channel_E = base_E.split('.')[2]
            
            #mtspec returns power spectra (square of spectra)
            stream = read(recordpath_N[0])
            tr = stream[0]
            data = tr.data
            delta = 1/tr.stats.sampling_rate
            print(delta)
            
            spec_amp, freq , jack, fstat, dof =  mtspec(data, delta = delta, time_bandwidth = 4, number_of_tapers=7, quadratic = True, statistics = True)
            
 			#find standard deviation
            sigmaN = (jack[:,1] - jack[:,0])/3.29
            
            #power spectra
            spec_array_N = np.array(spec_amp)
            freq_array_N = np.array(freq)
        
            stream = read(recordpath_E[0])
            tr = stream[0]
            data = tr.data

            spec_amp, freq , jack, fstat, dof =  mtspec(data, delta = delta, time_bandwidth = 4, number_of_tapers=7, quadratic = True, statistics = True)
            sigmaE = (jack[:,1] - jack[:,0])/3.29
            
            #power spectra
            spec_array_E = np.array(spec_amp)
            freq_array_E = np.array(freq)
            
            #if evenly sampled
            if(len(spec_array_E)==len(spec_array_N)):
                #here we bin into evenly spaced bins with frequency
                #spectra is power spectra so add the two components
                data_NE_2 = spec_array_E + spec_array_N
                #now data is NE power spectra
                #take the square root for normal velocity spectra
                data_NE = np.sqrt(data_NE_2)

                sigma = np.sqrt((spec_array_N/data_NE**2.)*sigmaN + ((spec_array_N/data_NE**2.)*sigmaN))
                
                #0.1-end
                bins, binned_data = bin_spec(data_NE[6:-1], freq[6:-1], num_bins = 75)
                bins_sig, binned_sig = bin_max_err(sigma[6:-1], freq[6:-1], num_bins = 75)

                #make sure that all spec is a number
                if (np.isnan(binned_data).any() == False):
                ##write to file
                    outfile = open(outpath + '/' +  event + '/' + network + '_' + station + '_' + 'HHNE' + '__' + event + '.out', 'w')
                    data = np.array([bins, binned_data, binned_sig])
#                    print(outpath + event)
                    data = data.T
                    outfile.write('#bins \t \t vel_spec_NE_m \t binned_sig \n')
                    np.savetxt(outfile, data, fmt=['%E', '%E', '%E'], delimiter='\t')
                    outfile.close()
    t2 = time.time()
    print('time for event: (s)', (t2-t1))
        
