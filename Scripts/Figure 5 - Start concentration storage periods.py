# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 16:29:43 2019

@author: emiel.kruisdijk
"""

#Open modules

import numpy as np
from matplotlib import pyplot as plt
import scipy
from pylab import *
import pandas as pd
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
import datetime
import functools
import math

########################################################################
#STYLE

#Legend
mpl.rcParams['legend.numpoints']=1
mpl.rcParams['legend.handlelength']=1
#mpl.legend(loc='best', frameon=True)
mpl.rcParams['legend.frameon'] = True
mpl.rcParams['legend.loc'] = 'best'
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['axes.labelsize'] = 6
mpl.rcParams['xtick.labelsize'] = 6
mpl.rcParams['ytick.labelsize'] = 6
mpl.rcParams['font.size'] = 6

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
########################################################################

df = pd.read_csv('../Data/Chem_data.csv', header= 0, skipinitialspace=True, sep=';' )

# Change dates to datetime dates

df['Date']=df['Date']+' '+df['Time']

#Makes from strings datetime formats
df['Date'] =  pd.to_datetime(df['Date'], format='%d-%m-%Y %H:%M')

chems=['IC_NO3 (mg/l)']

df['IC_NO3 (mg/l)'].where(df['IC_NO3 (mg/l)'] > 0.1, 0, inplace=True)

df1=df.loc[df['Project']=='BST2']
df2=df.loc[df['Project']=='BST4']
df3=df.loc[df['Project']=='BST5']

wells =['PZ3', 'PZ6', 'PZ1', 'PZ4','PZ5','PZ2']
wells2 =['MW 1', 'MW 2', 'MW 3', 'MW 4','MW 5','MW 6']
travtime = [38.0, 33.2, 21.9, 11.6, 25.0, 44.8]

fig, ax1= plt.subplots(nrows=1, ncols=1, figsize=(2.5, 2.5)) #figsize = A3
for tt, i,o in zip(travtime, wells, wells2):
    name1= df1.loc[df1['Name'] == i] 
    name2= df2.loc[df2['Name'] == i] 
    name3= df3.loc[df3['Name'] == i] 
    
    if i =='PZ6':
        ax1.plot(tt,name1['IC_NO3 (mg/l)'].values[0],  'o', c='k', label='Storage period - Winter 2019')
        ax1.plot(tt,name2['IC_NO3 (mg/l)'].values[0],  'o', c='b', label='Storage period - Fall 2020')
        ax1.plot(tt, name3['IC_NO3 (mg/l)'].values[0],  'o', c='r', label='Storage period - Spring 2021')
        ax1.set_xlabel('Min. traveltime (hours)')
        ax1.set_ylabel('NO$_\mathrm{3}$ concentration (mg/L)')
    else:
        ax1.plot(tt,name1['IC_NO3 (mg/l)'].values[0],  'o', c='k', label='_nolegend_')
        ax1.plot(tt,name2['IC_NO3 (mg/l)'].values[0],  'o', c='b', label='_nolegend_')
        ax1.plot(tt, name3['IC_NO3 (mg/l)'].values[0],  'o', c='r', label='_nolegend_')
        
ax1.annotate('MW4', (10.1, 5))
ax1.annotate('MW3', (20.4, 5))
ax1.annotate('MW5', (23.5, 21))
ax1.annotate('MW2', (31.7, 4.5))
ax1.annotate('MW1', (36.5, 18))
ax1.annotate('MW6', (42.5, 2))

ax1.legend(bbox_to_anchor=(0.5, -0.2), loc='upper center', ncol=1)
plt.savefig('../Figures/Figure_5.png',bbox_inches='tight', dpi=300)
