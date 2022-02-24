#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 10:02:31 2022
compare stress drops

@author: emmadevin
"""
import pandas as pd

path1 = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered/stress_drops/stress_drops_prelim_qa.csv'
path2 = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered7/stress_drops/stress_drops_prelim_qa.csv'


d1 = pd.read_csv(path1)
d2 = pd.read_csv(path2)

ev1 = d1['event']
sd1 = d1['stress drop']

ev2 = d2['event']
sd2 = d2['stress drop']

df1 = pd.DataFrame()
df2 = pd.DataFrame()

df1['Event'] = ev1
df1['Stress Drop Log Avg'] = sd1

df2['Event'] = ev2
df2['Stress Drop normal Avg'] = sd2


df3 = df1.merge(df2, how='inner', on='Event')