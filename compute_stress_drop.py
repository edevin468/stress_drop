#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02 Feb 2022

Estimating stress drop using Generalized Inversion Technique
Code written by Alexis Klimasewski and Emma Devin

INPUT: working directory for dataset (string), file type (string)
       working directory must contain the following: 

            ->dataset_name
                ->processed
                    ->ci38446071 (event id)
                        ->CI.APL.0.01_N.out (waveform files, in *.out or *.mseed format)
                ->RC_phase_beta (contains phase files for each event)
                    ->38446071.phase
                ->event_files (contains *.txt files for each event containing station metadata)

OUTPUT: directories for each step, plots where applicable, stress drops and other parameters (*.csv)

REQUIREMENTS: (in addition to usual python libraries)
    mtspec  
    obspy
    spec_func
    dread

"""
from stress_drop_fns import *

# define working directory name
working_dir = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered'

# define file type ('OUT' or 'MSEED')
file_type = 'OUT'

print('*******************************************')
print('Starting SETUP: compiling list of station information')
setup(working_dir)
print('SETUP COMPLETE.')
print('*******************************************')



print('*******************************************')
print('Starting STEP 1 of 5: computing velocity spectra and saving to record_spectra directory')
step1_compute_spectra(working_dir, file_type)
print('Step 1 of 5 COMPLETE.')
print('*******************************************')



print('*******************************************')
print('Starting STEP 2 of 5: inverting for event and site spectra and saving to Andrews_inversion directory')
step2_secondo_meters(working_dir)
print('Step 2 of 5 COMPLETE.')
print('*******************************************')



print('*******************************************')
print('Starting STEP 3 of 5: finding constraint event and Brune spectra and saving to constraint and Brune_spectra directories')
step3_findBrune_trap(working_dir)
print('Step 3 of 5 COMPLETE.')
print('*******************************************')



print('*******************************************')
print('Starting STEP 4 of 5: applying constraint function to all spectra and saving to Andrews_inversion_constrained directory')
step4_secondo_constraint(working_dir)
print('Step 5 of 5 COMPLETE.')
print('*******************************************')



print('*******************************************')
print('Starting STEP 5 of 5: finding best fit parameters and and saving to stress_drop directory')
step5_fitBrune(working_dir)
print('Step 5 of 5 COMPLETE.')
print('*******************************************')


print('Stress drop estimates successfully computed.')

