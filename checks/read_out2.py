#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 11:11:43 2021

@author: emmadevin
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob





 
    
# path = glob.glob('/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/record_spectra/ci38548295/*.out')
path = glob.glob('/Users/emmadevin/Work/USGS_2021/Data/Prelim+/record_spectra/38548295/*.out')
# path = glob.glob('/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/Andrews_inversion/Events/*.out')



for i in range(len(path)):
    data = np.genfromtxt(path[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1,2)) #only read in first two cols

    freq = data.T[0]
    spectra = data.T[1]
    spec_err = data.T[2]
    
    spec_err_log = np.abs(2*(data[:,2])/(data[:,1]*np.log(10)))
    
    
    
    fig = plt.figure(figsize = (4,4))
    plt.style.use('classic')
    fig.patch.set_facecolor('white')
    plt.plot(freq, spectra, lw=1,c= 'k')
    plt.errorbar(freq, spectra, yerr = spec_err, c='k')
    plt.title('', loc='left')
    # plt.xlim(0.04, 50)
    # plt.ylim(10**-8,10**2)
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('velocity amplitude (m)')
    
    
# path = glob.glob('/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered/Andrews_inversion/Events/*.out')

# for i in range(len(path)):
#     data = np.genfromtxt(path[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1,2)) #only read in first two cols

#     freq = data.T[0]
#     spectra = data.T[1]
    
    
#     plt.plot(freq, spectra, c= 'blue')
#     plt.title('(b)', loc='left')





# x = np.linspace(100,120,10)
# y = x

# plt.plot(x,y, c='grey', label='station spectra')
# plt.plot(x,y, c='blue', label='event spectra')
# plt.legend(loc = 'lower center', fontsize = 11)