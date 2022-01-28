#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 11:45:14 2021

@author: emmadevin
"""

import glob
import os.path as path
import numpy as np
import obspy
from obspy import read    
import dread
import time
import random
import pandas as pd
import matplotlib.pyplot as plt
fig = plt.figure()
fig.patch.set_facecolor('white')
plt.style.use('classic')

#read in the cut and corrected spectra = records
#records are in m/s
#should be binned frequencies and amplitudes
#make list of records and the corresponding events and stations

# working directory and outfil path for inversion
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered'
outfile_path = working_dir + '/Andrews_inversion'

# df with station locations
# use station_locs.csv to include uncorrected stations
stations = pd.read_csv(working_dir + '/Station_info/station_counts.csv') 

# df with station counts
counts = pd.read_csv(working_dir + '/Station_info/station_counts.csv')

#list of record files
ev = glob.glob(working_dir + '/record_spectra/*/*')

stn_id_list = []

# discard all records from stations with less that 3 records
for record in ev:
    filename = (record.split('/')[-1])
    base = path.basename(record)
    stn = base.split('_')[1]
    ntwk = base.split('_')[0]
    
    stns = list(counts['station'])
    ntwks = counts['network']
    
    index = stns.index(stn)
    count_list = counts['count']
    count = count_list[index]
    
 
    if count < 3: 
        ev.remove(record)
        
    if count >= 3:
       stn_id = ntwk+stn
       if stn_id not in stn_id_list:
           stn_id_list.append(stn_id)
        
record_path = ev
print('Number of records: ', len(record_path))
   

# get lists of station ids and station locations
stn_list = (stations['network']+stations['station']).tolist()
stn_lat = stations['latitude']
stn_lon = stations['longitude']


eventidlist = []
event_lat = []
event_lon = []
event_depth = []
record_freq = []
record_spec = []
record_std = []

t1 = time.time()

for i in range(len(record_path)):
    
    # read record filename and extract event and station identifiers
    record = (record_path[i].split('/')[-1])
    base = path.basename(record)
    eventid = base.split('_')[-1]
    eventid = eventid.split('.')[0]
    ntwk = base.split('_')[0]
    stn = base.split('_')[1]
    stn_id = ntwk + stn
    
    # print some updates to track progress
    print('Event: ', eventid)
    print('Station: ', stn)
    
    # read phase file fo get event locations
    phase_file = working_dir + '/RC_phase_beta/' + eventid + '.phase'
    phase = pd.read_csv(phase_file, sep = '\s+', index_col=0, nrows = 0).columns.tolist()
    
    # assign event coordinates and depth
    evlat = float(phase[4])
    evlon = float(phase[5])
    evdepth = float(phase[6])
    
    # assign station corrdinates and depth
    stlat = stn_lat[stn_list.index(stn_id)]
    stlon = stn_lon[stn_list.index(stn_id)]
    stdepth = 0
    
    #find distance between event and station
    dist =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
    
    #km to m
    dist = dist*1000.

    #read in spectra file
    data = np.genfromtxt(record_path[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1,2))  #only read in first two cols
    record_freq.append(data[:,0])
    
    # data is NE spectra; square for power spectra; correct for geometrical spreading 
    record_spec.append((data[:,1]*dist)**2.)
    
    # propagation here for errors going lin to log power
    std_prop = np.abs(2*(data[:,2])/(data[:,1]*np.log(10)))
    record_std.append(std_prop) #in log power here
    
    #if event info not part of lists yet add
    if eventid not in eventidlist:
        eventidlist.append(eventid)
        event_lat.append(evlat)
        event_lon.append(evlon)
        event_depth.append(evdepth)
        
t2 = time.time()

print('Time to read and distance correct all records: ', (t2-t1)/60.)
print('Number of records (frequency): ', len(record_freq))
print('Number of records (amplitude): ',len(record_spec))

freq_list = record_freq[0]
print(freq_list)
F_bins = len(freq_list)
print(F_bins)

rows = len(record_path) #testing first 10
print(rows)

index_matrix = [[0 for j in range(3)] for i in range(rows)]
stn_index_list = []

stn_test_list = []

for i in range(len(record_path)):
# for i in range(rows):
# for eventid in eventidlist:
    record = record_path[i].split('/')[-1]
    base = path.basename(record)
    
    print(base)
    
    network = base.split('_')[0]
    station = base.split('_')[1]
    stnid = network+station
    eventid = base.split('_')[-1]
    eventid = eventid.split('.')[0]
    
    print(network, station)
    print(eventid)
    
    stn_index = stn_list.index(stnid)
    if stn_index not in stn_index_list:
        stn_index_list.append(stn_index)
        
    if stnid not in stn_test_list:   
        stn_test_list.append(stnid)

    
    #make a tuple of record, event, station so indices can be assigned
    index_matrix[i] = [base, eventidlist.index(eventid), stn_id_list.index(stnid)]

# print(eventidlist[0])
# print(stn_list)

I = len(eventidlist)#events
J = len(stn_id_list)#stations
K = len(record_path)#records
K = rows



print('Number of events: ', I, ' Number of stations: ', J)
print('Number of rows (records): ', K, ' Number of cols (events+stations): ', I+J)

#make the G matrix of 1s and 0s and R matrix of records
G1 = np.zeros((K,I))
G2 = np.zeros((K,J))


for k in range(K):#for all records
    print(k)
    print(index_matrix[k][1])
    print(index_matrix[k][2])
    
    G1[k][index_matrix[k][1]] = 1 #record row, eventid col
    G2[k][index_matrix[k][2]] = 1 #record row, station col


G = np.concatenate((G1,G2), axis = 1)

print(G)

R = np.zeros((K,F_bins))
cov = np.zeros((K,F_bins))

#populate R matrix with the log of the record spectra (power record spectra)
for k in range(K):#for all records
    #each row is a record, and col is a frequency band
    #set row equal to the that spectral array
    #here we take the log
    R[k:,] = np.log10(record_spec[k])
    cov[k:,] = record_std[k]
    
m1 = np.zeros((I+J, F_bins))
m_cov = np.zeros((I+J, F_bins))

#do the inversion for each freq
for f in range(F_bins):
    t1 = time.time()
    d = R[:,f]#record for given frequency col
    dT = d.T
    print('inverting for frequency: ', f, freq_list[f])
    G_inv = np.linalg.pinv(G)
    covd = np.diag(cov[:,f])

    covm = np.dot((np.dot(G_inv, covd)), G_inv.T)
    m1[:,f] = np.dot(G_inv,dT)
    m_cov[:,f]= covm.diagonal()
    t2 = time.time()
    print('time for inversion: (min) ', round((t2-t1)/60., 4))


print(m1.shape)
#now split m into an event matrix and a station matrix
event = m1[0:I,:] #take first I rows
station = m1[I:I+J+1,:]
event_cov = m_cov[0:I,:]
station_cov = m_cov[I:I+J+1,:]
print(event.shape, station.shape)
print(event_cov.shape, station_cov.shape)

for i in range(I):#for each event
    print(i)
    print(eventidlist[i])
    print(outfile_path + '/' + eventidlist[i] + '.out')
    #go from the log of the power spectra to the regular spectra in m
    amp = np.sqrt(np.power(10.0, event[i,:]))

    std = (np.sqrt(np.abs(event_cov[i,:])/2.)*((amp)*(np.log(10))))
    
    filename = outfile_path + '/Events/' + eventidlist[i] + '.out'
    outfile = open(filename, 'w+')
    
    print(path.exists(filename))

    out = (np.array([freq_list, amp, std])).T
    
    outfile.write('#freq_bins \t vel_spec_NE_m \t stdev_m \n')
    np.savetxt(outfile, out, fmt=['%E', '%E', '%E'], delimiter='\t')
    outfile.close()
    
    # df = pd.DataFrame(out)
    # plt.plot(df[0],df[1])
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel('frequency (Hz)')
    # plt.ylabel('spectra')
    
    

fig = plt.figure()
fig.patch.set_facecolor('white')
plt.style.use('classic')
print('*************************')
for i in range(J):#for each station
    
    
    print(stn_list[i])
   
    amp = np.sqrt(np.power(10.0, station[i,:]))

    std1 = np.sqrt((station_cov[i,:]))
    std = np.abs((std1/2.)*(amp)*(np.log(10)))
    outfile = open(outfile_path + '/Stations/' + stn_list[i] + '.out', 'w')
    out = (np.array([freq_list, amp, std])).T
 
    outfile.write('#freq_bins \t vel_spec_NE_m \t stdev_m \n')
    np.savetxt(outfile, out, fmt=['%E', '%E', '%E'], delimiter='\t')
    outfile.close()
    
    df = pd.DataFrame(out)
    
    # plt.plot(df[0],df[1])
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel('frequency (Hz)')
    # plt.ylabel('spectra')
    # plt.xlim(10**-3,10**2)
    


    