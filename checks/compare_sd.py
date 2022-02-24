#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 10:02:31 2022
compare stress drops

@author: emmadevin
"""
import pandas as pd

OG_path = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_filtered/stress_drops/stress_drops_prelim_filtered.csv'
QA_path = '/Users/emmadevin/Work/USGS_2021/Data/Prelim_qa_filtered4/stress_drops/stress_drops_prelim_qa.csv'


og = pd.read_csv(OG_path)
qa = pd.read_csv(QA_path)

og_ev = og['event']
og_sd = og['stress drop']

qa_ev = qa['event']
qa_sd = qa['stress drop']

df1 = pd.DataFrame()
df2 = pd.DataFrame()

df1['Event'] = og_ev
df1['Stress Drop w/o QA'] = og_sd

df2['Event'] = qa_ev
df2['Stress Drop w/ QA'] = qa_sd


df3 = df1.merge(df2, how='inner', on='Event')