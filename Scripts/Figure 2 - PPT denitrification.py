# -*- coding: utf-8 -*-
#NOTE: Dashed green line is in wrong position! Get it in right position before going on!
"""
Created on Thu Mar  5 15:44:32 2020

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
from scipy.stats.stats import pearsonr

########################################################################
#STYLE

#Legend
mpl.rcParams['legend.numpoints']=1
mpl.rcParams['legend.handlelength']=1
#mpl.legend(loc='best', frameon=True)
mpl.rcParams['legend.frameon'] = True
mpl.rcParams['legend.loc'] = 'best'
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['axes.labelsize'] = 9
mpl.rcParams['xtick.labelsize'] = 9
mpl.rcParams['ytick.labelsize'] = 9
mpl.rcParams['font.size'] = 9

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
########################################################################

def Haggerty1998_PPT(chemdf, i, o, ax, r, tracer): #chem, r= compound, i = well
        global f, popt, sigma_ab, r_squared, pnr
        
        print(o, tracer)
        
        chemdf['tstart']=(chemdf['time2']-chemdf.loc[chemdf['Sample Text']== i+'-INF4', 'time2'].values[0])/ np.timedelta64(1, 'h')         #Time elapsed since the last sample, assuming that last sample is taking at end infiltration
        chemdf['tstart']=chemdf['tstart']/24
        Tinj=300/2/60 #Tinj in days, assuming a flow of 2 l/min
            
        #Only use reactant values, if Tracer is more than 20% and reactant above 0.1 mg/l!!! Mean of 4 INF samples is taken and last SP. Ct-GW/(INF-GW)
        chemdf['norm tracer']= (chemdf[tracer]-chemdf[tracer].loc[chemdf['Sample Text'].str.contains('GW')].values[0])/(chemdf[tracer].loc[chemdf['Sample Text'].str.contains('INF')].mean()-chemdf[tracer].loc[chemdf['Sample Text'].str.contains('GW')].values[0])
        chemdf['norm reactant']= (chemdf[r]-chemdf[r].loc[chemdf['Sample Text'].str.contains('GW')].values[0])/(chemdf[r].loc[chemdf['Sample Text'].str.contains('INF')].mean()-chemdf[r].loc[chemdf['Sample Text'].str.contains('GW')].values[0])

#        chemdf['norm tracer']= (chemdf[tracer]-chemdf[tracer].loc[chemdf['Sample Tet'].str.contains('SP14')].values[0])/(chemdf[tracer].loc[chemdf['Sample Tet'].str.contains('INF')].mean()-chemdf[tracer].loc[chemdf['Sample Tet'].str.contains('SP14')].values[0])
        chemdf=chemdf[chemdf['norm tracer'] > 0.2]
        #chemdf=chemdf[chemdf[r] > 0.1]
        
        chemdf = chemdf.loc[chemdf['Sample Text'].str.contains('SP')]
        
        x_values=np.arange(0,30, 0.1)
               
        try:
            chemdf['ln(Cr/Ctr)']=np.log(chemdf['norm reactant']/chemdf['norm tracer'])
            #chemdf['ln(Cr/Ctr)']=np.log((chemdf[r]/chemdf[r].values[0])/chemdf['norm tracer'])
            chemdf = chemdf[chemdf['ln(Cr/Ctr)'].notna()]
            chemdf = chemdf[chemdf['ln(Cr/Ctr)']!=-inf]
            
            if len(chemdf)>3:
                ax[1].plot(chemdf['tstart'], chemdf['ln(Cr/Ctr)'], "o", markersize=4, c='k', label='ln(C$_\mathrm{NO3}$/C$_\mathrm{tr}$)')
                if (o=='MW 4') or (o=='MW 5') or (o=='MW 6'):
                    ax[1].set_xlabel('Time (days)')
                if (o=='MW 1') or (o=='MW 4'):
                    ax[1].set_ylabel('ln(C$_\mathrm{NO3}$/C$_\mathrm{tr}$)')
                ax[1].set_title(o)
                #ax[1].set_ylim(ymax=1.8)
                
                #First-order, based on last three samples
                
               
                def func(x, a, b):
                    return a*x+b
                
                import scipy.optimize as optimization
                import scipy.stats as stats
                
###############################################################################
#Divide in before 5 days and after 5 days #####################################
###############################################################################
                before=chemdf.loc[chemdf['tstart']< 6]
                after = chemdf.loc[chemdf['tstart']> 5] 
                
                
                
                f=optimization.curve_fit(func, before['tstart'], before['ln(Cr/Ctr)'])[0][0]
                popt, pcov =optimization.curve_fit(func, before['tstart'], before['ln(Cr/Ctr)'])
                #95% confidence intervals
                sigma_ab = np.sqrt(np.diagonal(pcov)) #https://stackoverflow.com/questions/39434402/how-to-get-confidence-intervals-from-curve-fit
                
                residuals = before['ln(Cr/Ctr)']- func(before['tstart'], *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((before['ln(Cr/Ctr)']-np.mean(before['ln(Cr/Ctr)']))**2)
                r_squared = 1 - (ss_res / ss_tot)
                
                ax[1].plot(x_values, func(x_values, *popt), '--', c='k', label='k$_\mathrm{1}$='+str('%1.1e' %(popt[0]*-1))+' d$^\mathrm{-1}$')#+'\nr$^\mathrm{2}$='+str('%1.2f' %r_squared)) 
                #ax[1].plot(chemdf['tstart'], func(chemdf['tstart'], -0.01, popt[1]), '--', c='b', label='k=0.01')
                #ax[1].plot(chemdf['tstart'], func(chemdf['tstart'], -0.1, popt[1]), '--', c='r', label='k=0.1')
                #ax[1].plot(chemdf['tstart'], func(chemdf['tstart'], -0.25, popt[1]), '--', c='c', label='k=0.25')
                #ax[1].plot(chemdf['tstart'], func(chemdf['tstart'], -1, 0), '--', c='c')
                
                bound_upper = func(before['tstart'], *(popt + sigma_ab))
                bound_lower = func(before['tstart'], *(popt - sigma_ab))
                # plotting the confidence intervals
                #ax[1].fill_between(before['tstart'], bound_lower, bound_upper, color = 'black', alpha = 0.15)
                
                ax[1].legend(loc=3, labelspacing=0.3) 
                #ax[1].set_xlim(0,220)
                
                x=before['tstart']+8 #add value to increase all numbers above 0
                y=before['ln(Cr/Ctr)']
                pnr = pearsonr(x, y)
                
################################ AFTER ########################################
                
                f=optimization.curve_fit(func, after['tstart'], after['ln(Cr/Ctr)'])[0][0]
                popt, pcov =optimization.curve_fit(func, after['tstart'], after['ln(Cr/Ctr)'])
                #95% confidence intervals
                sigma_ab = np.sqrt(np.diagonal(pcov)) #https://stackoverflow.com/questions/39434402/how-to-get-confidence-intervals-from-curve-fit
                
                residuals = after['ln(Cr/Ctr)']- func(after['tstart'], *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((after['ln(Cr/Ctr)']-np.mean(after['ln(Cr/Ctr)']))**2)
                r_squared = 1 - (ss_res / ss_tot)
                
                ax[1].plot(x_values, func(x_values, *popt), ':', c='k', label='k$_\mathrm{2}$='+str('%1.1e' %(popt[0]*-1))+' d$^\mathrm{-1}$')#+'\nr$^\mathrm{2}$='+str('%1.2f' %r_squared)) 
                #ax[1].plot(chemdf['tstart'], func(chemdf['tstart'], -0.01, popt[1]), '--', c='b', label='k=0.01')
                #ax[1].plot(chemdf['tstart'], func(chemdf['tstart'], -0.1, popt[1]), '--', c='r', label='k=0.1')
                #ax[1].plot(chemdf['tstart'], func(chemdf['tstart'], -0.25, popt[1]), '--', c='c', label='k=0.25')
                #ax[1].plot(chemdf['tstart'], func(chemdf['tstart'], -1, 0), '--', c='c')
                
                bound_upper = func(after['tstart'], *(popt + sigma_ab))
                bound_lower = func(after['tstart'], *(popt - sigma_ab))
                # plotting the confidence intervals
                #ax[1].fill_between(after['tstart'], bound_lower, bound_upper, color = 'black', alpha = 0.15)
                
                ax[1].legend(loc=3, labelspacing=0.3) 
                ax[1].set_xlim(0,10)
                
                x=after['tstart']+8 #add value to increase all numbers above 0
                y=after['ln(Cr/Ctr)']
                pnr = pearsonr(x, y)
                
            else:
                print('Not enough measurements ', str(r))
                popt=[np.nan, np.nan]
                sigma_ab=[np.nan, np.nan]
                r_squared=np.nan
                pnr=[np.nan, np.nan]
                
        except:
            print('Not enough measurements v2', str(r))
            popt=[np.nan, np.nan]
            sigma_ab=[np.nan, np.nan]
            r_squared=np.nan
            pnr=[np.nan, np.nan]
            
        return chemdf
###############################################################################

chem = pd.read_csv('../Data/Chem_data.csv', header= 0, skipinitialspace=True, sep=';')
chem = chem.loc[chem['Project'].str.contains('PPT')]
chem = chem.loc[chem['Type']== 'Sample']
chem['time2'] =  pd.to_datetime(chem['Date']+' '+chem['Time'], format='%d-%m-%Y %H:%M')  

#delete outliers:
chem = chem[chem['Sample Text'] != 'PZ1-SP8']  


param = ['IC_NO3 (mg/l)']

tracer='IC_Br (mg/l)'

wells =['PZ3', 'PZ6', 'PZ1', 'PZ4','PZ5','PZ2']
wells2 =['MW 1', 'MW 2', 'MW 3', 'MW 4','MW 5','MW 6']

df=pd.DataFrame(columns= wells2)

for r in param:
    fig, axes= plt.subplots(nrows=2, ncols=3, figsize=(8, 4)) #figsize = A3
    for i,o, ax in zip(wells, wells2, enumerate(axes.flatten())):              
        if o != 'MW 5':
            tracer='IC_Br (mg/l)'
        else:
            tracer='IC_Cl (mg/l)'
            
        chemdf = chem.loc[chem['Name']== i]
       
        chemdf = Haggerty1998_PPT(chemdf, i, o, ax, r, tracer)
        
        if o== 'MW 1':
            ax[1].set_ylim(-1.1,0.2)
            ax[1].set_xlim(0,9)
        elif o== 'MW 2':
            ax[1].set_ylim(-4.5,0.5)
            ax[1].set_xlim(0,12)
        elif o== 'MW 3':
            ax[1].set_ylim(-2,0.5)
            ax[1].set_xlim(0,9)
        elif o== 'MW 4':
            ax[1].set_ylim(-0.5,0.2)
            ax[1].set_xlim(0,9)
        elif o== 'MW 5':
            ax[1].set_ylim(-0.4,0.2)
            ax[1].set_xlim(0,7)
        elif o== 'MW 6':
            ax[1].set_ylim(-2,0.5)
            ax[1].set_xlim(0,10)

        df.loc['NO3', o]=popt[0]*-1
        df.loc['NO3', '95% '+o]=sigma_ab[0]
        df.loc['NO3', 'r2 '+o]=r_squared
        df.loc['NO3', 'Rxy '+o]=pnr[0]
        df.loc['NO3', 'rho '+o]=pnr[1]
        
    plt.subplots_adjust(hspace=0.5, wspace=0.3)

    plt.savefig('../Figures/Figure_2.png',bbox_inches='tight', dpi=300)    

# Uncomment to get all statistical data in csv file.       
#df.to_csv('../Figures/Figure_2_raw_data.csv')




    
    