#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 09:36:08 2021

@author: emmadevin
"""

import numpy as np
import glob
import os.path as path
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, fit_report

#--------------------------------------------#
# BEGIN FUNCTION DEFINITIONS                 #
#--------------------------------------------#

def residual(params, freq, spec, sigma):
    beta = params['beta'].value
    U = params['U'].value
    rho = params['rho'].value
    fc = params['fc'].value
    M0 = params['M0'].value
    
    omega0 = (M0*U)/(4.*rho*np.pi*(beta**(3.0)))
    Brune_spec = (2.*np.pi*(freq)*omega0)/(1.+((1./fc)*freq)**2.)
    
    
    resids = Brune_spec - spec
    weighted = np.sqrt(resids ** 2 / sigma ** 2)
    
    return weighted

#--------------------------------------------#
# END FUNCTION DEFINITIONS                   #
#--------------------------------------------#
# fig = plt.figure(figsize = (6,6))
# plt.style.use('classic')
# fig.patch.set_facecolor('white')

working_dir =  '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered'
outfile_path = working_dir + '/stress_drops'

# define list of filepaths for all constrained event spectra
event_list = glob.glob(working_dir + '/Andrews_inversion_constrained/3*.out')

# define fixed Brune model parameters and add to params object for inversion
beta = 3500. #3500m/s
stressdrop = 5e6 #5e6 pascals
U = 0.63#0.63
rho = 2750. #2750kg/m^3

params = Parameters()
params.add('beta', value=3500, vary=False)
params.add('stressdrop', value=5e6, vary=False)
params.add('U', value=0.63, vary=False)
params.add('rho', value=2750, vary=False)

letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']

# obtain event id, magnitude, moment, Brune model corner frequency
id_list = []
mag_list = []
moment_list = []
Brune_fc_list = []
M_list = []

for event in event_list:
    
    # obtain event ids
    filename = (event.split('/')[-1])
    event_id = filename.split('.')[0]
    id_list.append(event_id)
    
    # obtain event magnitudes
    phase_file = working_dir + '/RC_phase_beta/' + event_id + '.phase'
    phase = pd.read_csv(phase_file, sep = '\s+', index_col=0, nrows = 0).columns.tolist()
    mag = float(phase[7])
    mag_list.append(mag)

    # if less than 3, convert local magnitude to moment magnitude
    if mag < 3.0:
        M = 0.884 + 0.754*mag  
    else:
        M = mag
        
    M_list.append(M)
    
    # moment from the moment magnitude
    M0 = 10.**((3./2.)*M + 9.1)
    moment_list.append(M0)
    
    # Brune corner frequency
    Brune_fc = beta*(stressdrop/(8.47*M0))**(1./3.)
    Brune_fc_list.append(Brune_fc)

# save all values computed above to dataframe
df = pd.DataFrame()

df['event'] = id_list
df['catalogue mag'] = mag_list
df['moment mag'] = M_list
df['catalogue moment'] = moment_list
df['Brune fc'] = Brune_fc_list 


# for each event preform nonlinear least squares to find fc and then compute stress drop
fc_list = []
sd_list = []
chisqr_list = []
M0_list = []
fc_std_list = []
M0_std_list = []
sd_std_list = []

for i in range(len(event_list)): 
    
    # parameters needed here
    Brune_fc = Brune_fc_list[i]
    params.add('M0', value = moment_list[i], vary=True)
    
    # read in data
    data = np.genfromtxt(event_list[i], dtype = float, comments = '#', delimiter = None)
    freq_untrimmed = np.array(data[:,0])
    spec_untrimmed = np.array(data[:,1])
    sigma_untrimmed = np.array(data[:,2])
    
    # Brune fc for this event
    Brune_fc = Brune_fc_list[i]
    params.add('fc', value=Brune_fc)
    
    # trim to desired frequency band
    upper= np.abs(freq_untrimmed - Brune_fc*5)
    upper_ind = upper.argmin()
    lower = np.abs(freq_untrimmed - Brune_fc/20)
    lower_ind = lower.argmin()
    spec = spec_untrimmed[lower_ind:upper_ind]
    freq = freq_untrimmed[lower_ind:upper_ind]
    sigma = sigma_untrimmed[lower_ind:upper_ind]
    
    # preform least squares inversion and assign to object
    out = minimize(residual, params, method='leastsq', args = (freq, spec, sigma))
   
    # make lists of the corner frequency best fit and chi^2 parameters
    fc = out.params['fc'].value
    fc_list.append(fc)
    
    fc_std = out.params['fc'].stderr
    fc_std_list.append(fc_std)
    
    M0 = out.params['M0'].value
    M0_list.append(M0)
    
    M0_std = out.params['M0'].stderr
    M0_std_list.append(M0_std)
    
    chisqr = out.chisqr
    chisqr_list.append(chisqr)
    
    # compute stress drop using equation from Shearer et al.
    sd = M0/(0.42*beta/fc)**3
    sd_list.append(sd/10**6)
    
    dsd_dM0 = fc**3/(0.42*beta)**3
    dsd_dfc = 3*M0*fc**2/(0.42*beta)**3
    
    sd_std = np.sqrt(dsd_dM0**2*M0_std**2 + dsd_dfc**2*fc_std**2)
    sd_std_list.append(sd_std/10**6)
    
    # make plot for each event to visually compare best fit curvlitudee with data
    beta = out.params['beta'].value
    U = out.params['U'].value
    rho = out.params['rho'].value
    fc = out.params['fc'].value
    
    
    omega0 = (M0*U)/(4.*rho*np.pi*(beta**(3.0)))
    Brune_spec = (2.*np.pi*(freq)*omega0)/(1.+((1./fc)*freq)**2.)
    Brune_spec_guess = (2.*np.pi*(freq)*omega0)/(1.+((1./Brune_fc)*freq)**2.)
    
    fig = plt.figure(figsize = (6,4))
    plt.style.use('classic')
    fig.patch.set_facecolor('white')
    plt.plot(freq, spec, c = 'b', label='constrained event spectra')
    plt.plot(freq, Brune_spec, c ='r', label = r'best-fit Brune model')
    plt.plot(freq, Brune_spec_guess, c ='gray', ls= '--',label = r'initial Brune model')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(freq[0],freq[-1])
    plt.ylabel('velocity amplitude (m)')
    plt.xlabel('frequncy (Hz)')
    # plt.title('(' + letters[i] + ') event: '+ id_list[i]+', M'+str(mag_list[i]), loc = 'left', fontsize = 11)
    plt.legend(loc='lower center', ncol=1)
    # plt.grid()
    
    plt.savefig(working_dir + '/fc_fitting_plots/' + id_list[i] + '.png', bbox_inches='tight')
    
    print(fit_report(out.params))
    
# save values to dataframe
df['chi^2'] = chisqr_list
df['best-fit fc'] = fc_list
df['best-fit M0'] = M0_list
df['stress drop'] = sd_list
df['std M0'] = M0_std_list
df['std fc'] = fc_std_list
df['std sd']  = sd_std_list


# df.to_csv(outfile_path + '/stress_drops_prelim_filtered.csv')
    
    