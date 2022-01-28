#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 15:14:47 2017

@author: Alexis Klimasewski

make a constraint file and put in the box directory
then constrain each event and station

if constraint is an event, divide event and multipy station
if constraint is station, multiply event and divide station

rename the outfile for the event or station

all in m/s
"""

import numpy as np
import glob
import os.path as path
import matplotlib.pyplot as plt


working_dir =  '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered'

#fill in constraint event
secondo_dir = '/constraint'
#fill in constraint event name here
constraint ='38548295'
constraint_file =  working_dir + secondo_dir + '/constraint_' + constraint + '.out'
brune_file = working_dir + '/Brune_spectra/' + constraint + '.out'
outfile_path = working_dir + '/Andrews_inversion_constrained'

con = np.genfromtxt(constraint_file)
cf_spec = con.T[1] 
freq_list = con.T[0]

brun = np.genfromtxt(brune_file)
brune_spec = brun.T[1]

ev_file = working_dir + '/Andrews_inversion/Events' + '/38548295.out'
ev = np.genfromtxt(ev_file)



fig = plt.figure(figsize = (5,4))
plt.style.use('classic')
fig.patch.set_facecolor('white')
plt.plot(con.T[0], con.T[1], c='green', lw = 2,label = 'constraint \nfunction')
plt.plot(ev.T[0],ev.T[1], c = 'blue',lw = 2,label = 'constraint \nevent spectra')
plt.plot(brun.T[0],brun.T[1], c='skyblue',lw = 4, label = 'constraint event \nbrune spectra')
plt.plot(ev.T[0], ev.T[1]/con.T[1], c = 'k',lw = 2, ls = '--',label = 'constraint \napplied to \nconstraint event')
plt.title('Event ID: 38548295, M4.9', loc='left')
plt.xlim(0.964784333/20,0.964784333*5)
plt.ylim(10**-1,10**2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('frequency (Hz)')
plt.ylabel('velocity amplitude (m)')
plt.legend(fontsize = 8,bbox_to_anchor=(1.2, -0.2),ncol= 4)
# plt.grid()


# secondo_ev =  glob.glob(working_dir + '/Andrews_inversion/Events' + '/*.out')
# secondo_stn = glob.glob(working_dir + '/Andrews_inversion/Stations' + '/*.out')

# ##not in log space anymore
# for i in range(len(secondo_ev)):#for each event
#     #make each row into an array
#     event = np.genfromtxt(secondo_ev[i])
#     eventid = path.basename(secondo_ev[i]).split('.')[0]
#     amp = event.T[1]/cf_spec
#     std = event.T[2]/cf_spec
    
#     if len(amp) == 0: 
#         print('no data event ' + eventid)

#     outfile = open(outfile_path + '/' + eventid + '.out', 'w')
#     out = (np.array([freq_list, amp, std])).T
#     outfile.write('#freq_bins \t vel_spec_NE_m \t stdev_m \n')
#     np.savetxt(outfile, out, fmt=['%E', '%E', '%E'], delimiter='\t')
#     outfile.close()
    
    
# for i in range(len(secondo_stn)):#for each station
#     #make each row into an array
#     stn = np.genfromtxt(secondo_stn[i])
#     stnid = path.basename(secondo_stn[i]).split('.')[0]
#     amp = stn.T[1]*cf_spec
#     std = stn.T[2]*cf_spec
    
    
#     outfile = open(outfile_path + '/' + stnid + '.out', 'w')
#     out = (np.array([freq_list, amp, std])).T
#     outfile.write('#freq_bins \t vel_spec_NE_m \t stdev_m \n')
#     np.savetxt(outfile, out, fmt=['%E', '%E', '%E'], delimiter='\t')

#     outfile.close()
    