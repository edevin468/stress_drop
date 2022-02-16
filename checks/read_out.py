#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 16:14:59 2022

@author: emmadevin
"""

import numpy as np
import matplotlib.pyplot as plt

# filepath = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/processed/ci38548295/CI.APL_0.01_Z.out'
# filepath = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/processed/ci38451079/CI.LRL_0.01_E.out'  
filepath = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/record_spectra/ci38451079/CI_CWC_NE__ci38451079.out'  


data = np.genfromtxt(filepath, dtype = float, comments = '#', delimiter = None) 

freq = data.T[0]
spectra = data.T[1]
spec_err = data.T[2]


fig = plt.figure(figsize = (4,4))
plt.style.use('classic')
fig.patch.set_facecolor('white')
plt.plot(freq, spectra, lw=3,c= 'k')
plt.errorbar(freq, spectra,yerr = spec_err , c='k')
plt.title('', loc='left')
plt.xlim(0.04, 50)
# plt.ylim(10**-4,10**2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('frequency (Hz)')
plt.ylabel('velocity amplitude (m)')