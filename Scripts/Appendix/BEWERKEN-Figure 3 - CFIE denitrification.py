# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 12:01:34 2021

@author: emiel.kruisdijk
"""
import numpy as np
from matplotlib import pyplot as plt
import scipy
from pylab import *
import pandas as pd
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
import datetime
import functools
import matplotlib.patches as patches
import sys
import math
from scipy import special
import scipy.optimize as optimization

########################################################################
#STYLE

#Legend
mpl.rcParams['legend.numpoints']=1
mpl.rcParams['legend.handlelength']=1
#mpl.legend(loc='best', frameon=True)
mpl.rcParams['legend.frameon'] = True
mpl.rcParams['legend.loc'] = 'best'
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['font.size'] = 8

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
########################################################################

def Mean_Inf(chem):
    name=chem.loc[chem['Project']=='BST1 - INF']
    
    chems2 = [s for s in name.columns.values if 'Compound' in s]
    chems = [s for s in chems2 if '-d' not in s]
    
    chems2=np.append(chems2, 'IC_Br (mg/l)')
    chems2=np.append(chems2, 'IC_Cl (mg/l)')
    chems2=np.append(chems2, 'IC_NO3 (mg/l)')
    
    average={}
    for i in chems2:
        b=name[i].mean()
        average[i] = b 


    maxy={}
    for i in chems2:
        b=name[i].quantile(.9)
        maxy[i] = b   
        
    miny={}
    for i in chems2:
        b=name[i].quantile(.1)
        miny[i] = b
        
    return average, maxy, miny

def data_processing(chemdf, i):
    chem = chemdf.loc[chemdf['Type']== 'Sample']
    chem = chem.loc[chem['Project'].str.contains('BST1')]
    chem['time2'] =  pd.to_datetime(chem['Date']+' '+chem['Time'], format='%d-%m-%Y %H:%M')
    chem = chem.loc[chem['Name']== i]
    chem=chem.sort_values(by=['time2'])
    chem = chem.loc[~chem['Sample Text'].str.contains('NGW')]

    return chem

wells =['PZ3', 'PZ6', 'PZ1', 'PZ4','PZ5','PZ2']
wells2 =['MW-1', 'MW-2', 'MW-3', 'MW-4','MW-5','MW-6']
depths =['11.4-12.4 m-b.s.l.','15.0-16.0 m-b.s.l.','18.3-19.3 m-b.s.l.','22.8-23.8 m-b.s.l.','25.9-26.9 m-b.s.l.','31.2-32.2 m-b.s.l.']

############### Total infiltrated volume ######################################
System_start=['01-11-2019 09:21', '02-11-2019 08:40', '03-11-2019 08:50', '04-11-2019 7:42', '05-11-2019 07:50', '06-11-2019 09:45']
System_end=['01-11-2019 12:58', '02-11-2019 19:10', '03-11-2019 17:05', '04-11-2019 18:10', '05-11-2019 17:20', '06-11-2019 16:45']
Volume_start_INF1=[31.73, 39.71, 137.60, 218.16, 319.69, 408.24]
Volume_start_INF2=[29.35, 40.89, 137.29, 214.43, 311.08, 395.29]
Volume_end_INF1=[39.59, 137.60, 218.16, 319.69, 408.24, 471.54]
Volume_end_INF2=[40.79, 139.27, 214.43, 311.08, 395.29, 454.83]

d={'System_start':System_start, 'System_end':System_end, 'Volume_start_INF1':Volume_start_INF1, 'Volume_start_INF2':Volume_start_INF2, 'Volume_end_INF1':Volume_end_INF1, 'Volume_end_INF2':Volume_end_INF2}

System=pd.DataFrame(data=d)

System['System_start'] =  pd.to_datetime(System['System_start'], format='%d-%m-%Y  %H:%M')
System['System_end'] =  pd.to_datetime(System['System_end'], format='%d-%m-%Y  %H:%M')


chemdf = pd.read_csv('../Data/Chem_data.csv', header= 0, skipinitialspace=True, sep=';')
average, maxy, miny=Mean_Inf(chemdf)

fig, axes= plt.subplots(nrows=2, ncols=3, figsize=(8, 4)) #figsize = A3
for i,o, ax in zip(wells, wells2, enumerate(axes.flatten())): 
    chem = data_processing(chemdf, i)
    
    ax[1].plot(chem['time2'], chem['IC_NO3 (mg/l)'], 'o', c='k', markersize=5, label='Observed NO$_\mathrm{3}$ (mg/l)', zorder=100)
    
    for index, row in System.iterrows():
        rect = patches.Rectangle((row['System_start'],-100000),(row['System_end']-row['System_start']),300000,linewidth=1,edgecolor='k',facecolor='k', alpha=0.1)
        ax[1].add_patch(rect)
        
    ax[1].set_xlim(xmax=pd.to_datetime('19-Nov-2019'))
    ax[1].set_ylim(0,13)
    ax[1].tick_params(axis='x', rotation=90)
    ax[1].axhline(y=average['IC_NO3 (mg/l)'], color='b', linestyle='--', label='Mean inj. TDW NO$_\mathrm{3}$')
    ax[1].axhline(y=miny['IC_NO3 (mg/l)'], color='b', linestyle='--', label='Min inj. TDW NO$_\mathrm{3}$')
    ax[1].axhline(y=maxy['IC_NO3 (mg/l)'], color='b', linestyle='--', label='Max inj. TDW NO$_\mathrm{3}$') #only for legend
    ax[1].axhline(y=1, color='w', linestyle='--', label=' ', zorder=0)
    ax[1].set_title(o)
    
    ax2=ax[1].twinx()
    ax2.plot(chem['time2'], chem['IC_Cl (mg/l)'], 'o', c='k', markersize=5, alpha=0.2, label='Observed Cl')
    
    if o == 'MW-6':
        df=pd.read_csv('../Data/PHREEQC/MW6.txt', header= 0, skipinitialspace=True, sep='\t')
        df['Day_td']=pd.to_timedelta(df['Day'], unit='h')
        df['days2']=df['Day_td']+pd.to_datetime('01-11-2019 00:00', format='%d-%m-%Y  %H:%M')
        ax2.plot(df['days2'], df['IC_Cl_(mg/l)'], '-', c='k', alpha=0.2)
        ax[1].plot(df['days2'], df['NO3_0.0'], '-', c='b')
        ax[1].plot(df['days2'], df['NO3_0.1'], '-', c='c')
        ax[1].plot(df['days2'], df['NO3_0.25'], '-', c='m')
        ax[1].plot(df['days2'], df['NO3_0.5'], '-', c='tab:orange')
        ax[1].plot(df['days2'], df['NO3_1.0'], '-', c='tab:purple')
    elif o == 'MW-5':
        df=pd.read_csv('../Data/PHREEQC/MW5.txt', header= 0, skipinitialspace=True, sep='\t')
        df['Day_td']=pd.to_timedelta(df['Day'], unit='h')
        df['days2']=df['Day_td']+pd.to_datetime('01-11-2019 00:00', format='%d-%m-%Y  %H:%M')
        ax2.plot(df['days2'], df['IC_Cl_(mg/l)'], '-', c='k', alpha=0.2)
        ax[1].plot(df['days2'], df['NO3_0.0'], '-', c='b')
        ax[1].plot(df['days2'], df['NO3_0.1'], '-', c='c')
        ax[1].plot(df['days2'], df['NO3_0.25'], '-', c='m')
        ax[1].plot(df['days2'], df['NO3_0.5'], '-', c='tab:orange')
        ax[1].plot(df['days2'], df['NO3_1.0'], '-', c='tab:purple')
    elif o == 'MW-4':
        df=pd.read_csv('../Data/PHREEQC/MW4.txt', header= 0, skipinitialspace=True, sep='\t')
        df['Day_td']=pd.to_timedelta(df['Day'], unit='h')
        df['days2']=df['Day_td']+pd.to_datetime('01-11-2019 00:00', format='%d-%m-%Y  %H:%M')
        ax2.plot(df['days2'], df['IC_Cl_(mg/l)'], '-', c='k', alpha=0.2, label='Simulated Cl')
        ax[1].plot(df['days2'], df['NO3_0.0'], '-', c='b', label='k= 0.0 d$^\mathrm{-1}$') #NO$_\mathrm{3}$ (mg/l)
        ax[1].plot(df['days2'], df['NO3_0.1'], '-', c='c', label='k= 0.1 d$^\mathrm{-1}$')
        ax[1].plot(df['days2'], df['NO3_0.25'], '-', c='m', label='k= 0.25 d$^\mathrm{-1}$')
        ax[1].plot(df['days2'], df['NO3_0.5'], '-', c='tab:orange', label='k= 0.5 d$^\mathrm{-1}$')
        ax[1].plot(df['days2'], df['NO3_1.0'], '-', c='tab:purple', label='k= 1.0 d$^\mathrm{-1}$')
        ax2.axhline(y=1, color='w', linestyle='--', label=' ', alpha=0)
        ax2.axhline(y=1, color='w', linestyle='--', label=' ', alpha=0)
        ax2.axhline(y=1, color='w', linestyle='--', label=' ', alpha=0)
        handles1, labels1 = ax[1].get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        
        handles=handles1+handles2
        labels=labels1+labels2
        ax[1].legend(handles, labels, loc='upper center', bbox_to_anchor=(2.0, -0.8), ncol=3)
        
    elif o == 'MW-3':
        ax[1].tick_params(axis='x',which='both',bottom='on', top='off', labelbottom='off')
        ax2.tick_params(axis='x',which='both',bottom='on', top='off', labelbottom='off')
        df=pd.read_csv('../Data/PHREEQC/MW3.txt', header= 0, skipinitialspace=True, sep='\t')
        df['Day_td']=pd.to_timedelta(df['Day'], unit='h')
        df['days2']=df['Day_td']+pd.to_datetime('01-11-2019 00:00', format='%d-%m-%Y  %H:%M')
        ax2.plot(df['days2'], df['IC_Cl_(mg/l)'], '-', c='k', alpha=0.2)
        ax[1].plot(df['days2'], df['NO3_0.0'], '-', c='b')
        ax[1].plot(df['days2'], df['NO3_0.1'], '-', c='c')
        ax[1].plot(df['days2'], df['NO3_0.25'], '-', c='m')
        ax[1].plot(df['days2'], df['NO3_0.5'], '-', c='tab:orange')
        ax[1].plot(df['days2'], df['NO3_1.0'], '-', c='tab:purple')
    elif o == 'MW-2':
        ax[1].tick_params(axis='x',which='both',bottom='on', top='off', labelbottom='off')
        ax2.tick_params(axis='x',which='both',bottom='on', top='off', labelbottom='off')
        df=pd.read_csv('../Data/PHREEQC/MW2.txt', header= 0, skipinitialspace=True, sep='\t')
        df['Day_td']=pd.to_timedelta(df['Day'], unit='h')
        df['days2']=df['Day_td']+pd.to_datetime('01-11-2019 00:00', format='%d-%m-%Y  %H:%M')
        ax2.plot(df['days2'], df['IC_Cl_(mg/l)'], '-', c='k', alpha=0.2)
        ax[1].plot(df['days2'], df['NO3_0.0'], '-', c='b')
        ax[1].plot(df['days2'], df['NO3_0.1'], '-', c='c')
        ax[1].plot(df['days2'], df['NO3_0.25'], '-', c='m')
        ax[1].plot(df['days2'], df['NO3_0.5'], '-', c='tab:orange')
        ax[1].plot(df['days2'], df['NO3_1.0'], '-', c='tab:purple')
    elif o == 'MW-1':
        ax[1].tick_params(axis='x',which='both',bottom='on', top='off', labelbottom='off')
        ax2.tick_params(axis='x',which='both',bottom='on', top='off', labelbottom='off')
        df=pd.read_csv('../Data/PHREEQC/MW1.txt', header= 0, skipinitialspace=True, sep='\t')
        df['Day_td']=pd.to_timedelta(df['Day'], unit='h')
        df['days2']=df['Day_td']+pd.to_datetime('01-11-2019 00:00', format='%d-%m-%Y  %H:%M')
        ax2.plot(df['days2'], df['IC_Cl_(mg/l)'], '-', c='k', alpha=0.2)
        ax[1].plot(df['days2'], df['NO3_0.0'], '-', c='b')
        ax[1].plot(df['days2'], df['NO3_0.1'], '-', c='c')
        ax[1].plot(df['days2'], df['NO3_0.25'], '-', c='m')
        ax[1].plot(df['days2'], df['NO3_0.5'], '-', c='tab:orange')
        ax[1].plot(df['days2'], df['NO3_1.0'], '-', c='tab:purple')

fig.text(0.05,0.7, 'NO$_\mathrm{3}$ concentrations (mg/L)', rotation='vertical', fontsize=8)
fig.text(0.98,0.7, 'Cl concentrations (mg/L)', rotation='vertical', fontsize=8)
plt.subplots_adjust(hspace=0.3, wspace=0.5)
plt.savefig('../Figures/Figure_3.png',bbox_inches='tight', dpi=300)    
