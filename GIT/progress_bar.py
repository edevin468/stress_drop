#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 09:37:54 2021

progress bar

@author: emmadevin
"""

import sys

class progress_bar:
    def __init__(self, l=0):
        self.length = l
    

    
    def get_progress(self, i):
        self.step = i+1
        prog = self.step/self.length
        
        stars = int(prog*50)
        
        spaces = int(50-(prog*50))
        
        
        if stars + spaces != 50: stars+=1
        
        
        fill_list = ['#']*stars
        
        fill = ''
        for e in fill_list:
            fill+=e
        
        
        space = ''
        space_list = [' ']*spaces
        for e in space_list:
            space+=e
            
        
        
        bar = fill + space
    
            
        sys.stdout.write('\rProcessing: |' + bar + '| (' + f'{prog*100:.2f}' +'%)')
        sys.stdout.flush()
        
        
    

