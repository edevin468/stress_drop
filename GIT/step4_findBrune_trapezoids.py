#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 09:53:15 2017

@author: Alexis Klimasewski

compute “B” for a M3 earthquake, and then for every ~M3 earthquake
in your dataset (call it “A”), find A/B for every frequency range you inverted for.  
Then sum up A/B over all frequency ranges and find which earthquake has the minimum value for that
because then that would suggest that earthquake is closest to a brune spectrum
"""

import glob
import numpy as np
import obspy
import matplotlib.pyplot as plt
import pandas as pd

#properties for unit
beta = 3500. #3500m/s
stressdrop = 5e6 #5e6 pascals
U = 0.63#0.63
rho = 2750. #2750kg/m^3

working_dir =  '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered'
event_spectra_dir = working_dir + '/Andrews_inversion/Events/'
event_list = glob.glob(event_spectra_dir + '*.out')

writefile = 'yes'


    
#compute Brune spectra for all of the events in directory
cf_list = []
cf2_list = []
Brune_list = []
Brune_list_untrimmed = []
spec_list = []
ev_list = []
mag_list = []
freq_list = []
freq_untrimmed_list = []
id_list = []

spec_demean_list = []
Brune_demean_list = []

spec_demean_list_log = []
Brune_demean_list_log = []

event_spectra = []

for event in event_list:
    
    filename = (event.split('/')[-1])
    event_id = filename.split('.')[0]
    
    # obtain event magnitudes
    phase_file = working_dir + '/RC_phase_beta/' + event_id + '.phase'
    phase = pd.read_csv(phase_file, sep = '\s+', index_col=0, nrows = 0).columns.tolist()
    
    mag = float(phase[7])
    
    m_l = 1
    m_u = 5
    
    if mag <= m_u and mag >= m_l:
        keep = True
    else:
        keep = False
        
    if keep == True:
        event_spectra.append(event)
    else:
        continue
    

    
for event in event_spectra:
    
    filename = (event.split('/')[-1])
    event_id = filename.split('.')[0]
    print(event_id)
    id_list.append(event_id)
    
    # obtain event magnitudes
    phase_file = working_dir + '/RC_phase_beta/' + event_id + '.phase'
    phase = pd.read_csv(phase_file, sep = '\s+', index_col=0, nrows = 0).columns.tolist()
    
    mag = float(phase[7])
    print(mag)
    mag_list.append(mag)
    
    
    data = np.genfromtxt(event, dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
    freq_untrimmed = np.array(data[:,0])
    spec_untrimmed = np.array(data[:,1])  # these record spectra are in m

    
    # if less than 3, convert local magnitude to moment magnitude
    if mag < 3.0:
        M = 0.884 + 0.754*mag  # 0.884 + 0.667*ml, 754
    else:
        M = mag
        
    # compute Brune in SI units
    # moment from the moment magnitude
    M0 = 10.**((3./2.)*M + 9.1)
    
    #corner frequency
    fc = beta*(stressdrop/(8.47*M0))**(1./3.)
    omega0 = (M0*U)/(4.*rho*np.pi*(beta**(3.0)))
    
    freq_untrimmed_list.append(freq_untrimmed)
   
    #brune spectra over all frequencies
    Brune_untrimmed = (2.*np.pi*(freq_untrimmed)*omega0)/(1.+((1./fc)*freq_untrimmed)**2.)
    
    upper= np.abs(freq_untrimmed - fc*5)
    upper_ind = upper.argmin()
    
    lower = np.abs(freq_untrimmed - fc/20)
    lower_ind = lower.argmin()
    
    spec = spec_untrimmed[lower_ind:upper_ind]
    freq = freq_untrimmed[lower_ind:upper_ind]
    
    #brune spectra over all frequencies
    Brune = (2.*np.pi*(freq)*omega0)/(1.+((1./fc)*freq)**2.)


    Brune_list.append(Brune)
    Brune_list_untrimmed.append(Brune_untrimmed)
    spec_list.append(spec)
    freq_list.append(freq)
    
    #stay in meters
    shift1 = np.mean(Brune[0:74])
    shift2 = np.mean(spec[0:74])
    
    # cf_list.append(np.log10(spec_untrimmed/shift2)-np.log10(Brune_untrimmed/shift1))
    cf_list.append(np.log10(spec_untrimmed)-np.log10(Brune_untrimmed))
    
cfarray = np.array(cf_list)
    

# find difference between areas under spectral curve and Brune curve using area between curves
areas = []
print('***********')
for i in range(len(spec_list)):
    spectra = spec_list[i]
    brune = Brune_list[i]
    freq = freq_list[i]
    area = 0
    
    # normalize data and brune by avg amplitude
    avg_spec = np.mean(spectra)
    avg_brune = np.mean(brune)
    
    # spectra = 10**(np.log10(spectra)-avg_spec)
    # brune = 10**(np.log10(brune)-avg_brune)
    
    spectra = spectra/avg_spec
    brune = brune/avg_brune
    
    print('M',mag_list[i])
    print(np.mean(brune))
    print(np.mean(spectra))
    
    # find difference in amplitude and remove it to only comepare shapes
    # amp_diff = np.mean(spectra) - np.mean(brune)
    # print(amp_diff)
    
    # if amp_diff > 0:
    #     brune = 10**(np.log10(brune)+ amp_diff)
    # elif amp_diff < 0:
    #     brune = 10**(np.log10(brune)- amp_diff)
    # else: 
    #     brune = brune
    
    
    for j in range(len(spectra)-1): 
        spec_trap = 1/2*(freq[j+1]-freq[j])*(spectra[j+1]+spectra[j])
        Brune_trap =  1/2*(freq[j+1]-freq[j])*(brune[j+1]+brune[j])
        
        area += abs(spec_trap-Brune_trap)
        
    areas.append(area)
    
    # print('_____________________')
    # print('M',mag_list[i])
    # print(area)
    # print(np.mean(brune))
    # print(np.mean(spectra))
    
    
    fig = plt.figure(figsize = (8,6))
    plt.style.use('classic')
    fig.patch.set_facecolor('white')
    
    plt.ylabel('velocity amplitude (m)', fontsize = 12)
    # plt.xlim(0.01,100)
    plt.loglog(freq , spectra, color = 'green', label = 'event spectra')
    plt.grid()
    plt.loglog(freq, brune, color = 'blue', label = 'Brune spectra')
    plt.legend(loc = 'lower left', fontsize = 12)
    plt.xlabel('frequency (Hz)', fontsize = 12)
    # plt.title('event ID: 38475431 \nmagnitude : 4.15)
    plt.title('event id:'+id_list[i]+', magnitude: M'+str(mag_list[i]), fontsize = 12, loc='left')
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    
    

diff_list = list(areas)
ind = diff_list.index(min(diff_list))



print(event_spectra[ind])
print(min(diff_list))        

print(mag_list[ind])
fig = plt.figure(figsize = (8,6))
plt.style.use('classic')
fig.patch.set_facecolor('white')


# plt.ylabel('velocity amplitude (m)', fontsize = 12)
# # plt.xlim(0.01,100)
# plt.loglog(freq_list[ind] , spec_list[ind], color = 'green', label = 'event spectra')
# plt.grid()
# plt.loglog(freq_list[ind], Brune_list[ind], color = 'blue', label = 'Brune spectra')
# plt.legend(loc = 'lower left', fontsize = 12)
# plt.xlabel('frequency (Hz)', fontsize = 12)
# # plt.title('event ID: 38475431 \nmagnitude : 4.15)
# plt.title('event id:'+id_list[ind]+', magnitude: M'+str(mag_list[ind]), fontsize = 12, loc='left')
# plt.tick_params(axis='both', which='major', labelsize=12)
# plt.tick_params(axis='both', which='both', length = 5, width = 1)
# plt.text(0.7, .1, 'Median log(diff) 1-32.7 Hz (demeaned): ' + str(round(sum_list[ind],3)), fontsize = 16)

# plt.savefig('/Users/emmadevin/Work/USGS 2021/Figures/SCEC/Brune_example.pdf')
# plt.show()



#write the constraint file in linear space to agree with the event and station spectra
if writefile == 'yes':
    outfile = open(working_dir + '/constraint/constraint_' + id_list[ind] + '.out', 'w')
    out = np.array([freq_untrimmed_list[ind],10.**(cf_list[ind])]).T
    outfile.write('#freq_bins \t cf_m \n')
    np.savetxt(outfile, out, fmt='%E', delimiter='\t')
    outfile.close()
    
for ind in range(len(Brune_list_untrimmed)):
    outfile = open(working_dir + '/Brune_spectra/' + id_list[ind] + '.out', 'w')
    out = np.array([freq_untrimmed_list[ind],Brune_list_untrimmed[ind]]).T
    outfile.write('#freq_bins \t brune \n')
    np.savetxt(outfile, out, fmt='%E', delimiter='\t')
    outfile.close()

