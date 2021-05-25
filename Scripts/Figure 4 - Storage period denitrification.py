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

rates=pd.DataFrame()

for r in chems:
#    try:
        fig, axes= plt.subplots(nrows=2, ncols=3, figsize=(5, 2.5)) #figsize = A3
        plt.subplots_adjust(hspace=0.3, wspace=0.3)
        for i,o, ax in zip(wells, wells2, enumerate(axes.flatten())):
            #BST2
            name1 = df1.loc[df['Name']== i]
            name1=name1.sort_values('Date')
            name1['hours']=(name1['Date']-name1['Date'].values[0])/ np.timedelta64(1, 'h')
            
            name1 = name1[name1[r].notna()]
            name1 = name1[name1[r]!=-inf]
            name1['days']=name1['hours']/24
            #BST4
            name2 = df2.loc[df['Name']== i]
            name2=name2.sort_values('Date')
            name2['hours']=(name2['Date']-name2['Date'].values[0])/ np.timedelta64(1, 'h')
            
            name2 = name2[name2[r].notna()]
            name2 = name2[name2[r]!=-inf]
            name2['days']=name2['hours']/24
            #BST5
            name3 = df3.loc[df['Name']== i]
            name3=name3.sort_values('Date')
            name3['hours']=(name3['Date']-name3['Date'].values[0])/ np.timedelta64(1, 'h')
            
            name3 = name3[name3[r].notna()]
            name3 = name3[name3[r]!=-inf]
            name3['days']=name3['hours']/24
            if len(name1) > 3 and len(name2) > 3 and len(name3) > 3:
###############################################################################
#-----------------------plot observations-------------------------------------#
###############################################################################
                

                if 'Compound' in r:
                    for i in bieb:
                        if i in r:
                            ax[1].axhline(y=bieb[i][1], color='g', linestyle='--', alpha=0.5, label='_nolegend_')
                            ax[1].plot(name1['days'], name1[r], "o", markersize=3, c='k', label='_nolegend_') # multiple line can be done with \n
                            #ax[1].errorbar(name1['days'], name1[r], yerr=0.2*name1[r], fmt='o', c='k',markersize=1, linewidth=0.7, label='_nolegend_')   
                            
                            ax[1].axhline(y=bieb[i][1], color='g', linestyle='--', alpha=0.5, label='_nolegend_')
                            ax[1].plot(name2['days'], name2[r], "o", markersize=3, c='b', label='_nolegend_') # multiple line can be done with \n
                            #ax[1].errorbar(name2['days'], name2[r], yerr=0.2*name2[r], fmt='o', c='b',markersize=1, linewidth=0.7, label='_nolegend_')
                            
                            ax[1].axhline(y=bieb[i][1], color='g', linestyle='--', alpha=0.5, label='_nolegend_')
                            ax[1].plot(name3['days'], name3[r], "o", markersize=3, c='r', label='_nolegend_') # multiple line can be done with \n
                            #ax[1].errorbar(name2['days'], name2[r], yerr=0.2*name2[r], fmt='o', c='b',markersize=1, linewidth=0.7, label='_nolegend_')  

                else:
                    ax[1].plot(name1['days'], name1[r], "o", markersize=2, c='k', label='_nolegend_') # multiple line can be done with \n
                    ax[1].plot(name2['days'], name2[r], "o", markersize=2, c='b', label='_nolegend_') # multiple line can be done with \n
                    ax[1].plot(name3['days'], name3[r], "o", markersize=2, c='r', label='_nolegend_') # multiple line can be done with \n
                    
                ax[1].set_title(o)
                
                if (o == 'MW 4') or (o == 'MW 5') or (o == 'MW 6'):
                    ax[1].tick_params(axis='x', rotation=0)
                else:
                    ax[1].tick_params(axis='x',which='both',bottom='on', top='off', labelbottom='off')
###############################################################################                
#--------------------plot fit + uncertainties---------------------------------#
###############################################################################
                if name1[r].max()>0.1:
                    x0=0.01
                    def func(x, k, b):
                        return np.exp(-k*x)*b
                    import scipy.optimize as optimization
                    #BST2
                    #Initial guess
                    k0=0.001
                    A0=name1[r].values[0]
                    name1['days']=name1['hours']/24
                    popt, pcov =optimization.curve_fit(func, name1['days'], name1[r])
                    #95% confidence intervals
                    sigma_ab = np.sqrt(np.diagonal(pcov)) #https://stackoverflow.com/questions/39434402/how-to-get-confidence-intervals-from-curve-fit   
                    residuals = name1[r]- func(name1['days'], *popt)
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((name1[r]-np.mean(name1[r]))**2)
                    r_squared = 1 - (ss_res / ss_tot)
                    dt50=np.log(2)/popt[0]
                    dt50_95=np.log(2)/sigma_ab[0]
                    
                    timestep=np.arange(0,50,0.1)
                    ax[1].plot(timestep, func(timestep, *popt), '--', c='k', linewidth=0.7, label=str('%4.2f' %(popt[0])+'d$^\mathrm{-1}$'))#+', r$^\mathrm{2}$= '+('%1.2f' %(r_squared)))#+'    +-'+str('%.1E' %(sigma_ab[0]))) #'\nDT50='+str('%4.1f' %(np.log(2)/popt[0]))+'d') 
              
                    #bound_upper = func(name1['days'], (popt[0]+sigma_ab[0]), (popt[1]-sigma_ab[1]))
                    #bound_lower = func(name1['days'], (popt[0]-sigma_ab[0]), (popt[1]+sigma_ab[1]))
                    # plotting the confidence intervals
                    #ax[1].fill_between(name1['days'], bound_lower, bound_upper,
                    #                  color = 'black', alpha = 0.15)
                    
                    if 'Compound'in r:
                        for i in bieb:
                            if i in r:
                                rates.loc['BST2 '+bieb[i][2],'DT50 '+ o]=dt50
                                rates.loc['BST2 '+bieb[i][2],'DT50 95% '+ o]=dt50_95
                                rates.loc['BST2 '+bieb[i][2], 'k '+ o]=popt[0]
                                rates.loc['BST2 '+bieb[i][2], 'A0 '+ o]=popt[1]
                                rates.loc['BST2 '+bieb[i][2], 'k95% '+o]=sigma_ab[0] 
                                rates.loc['BST2 '+bieb[i][2], 'A095% '+o]=sigma_ab[1] 
                                rates.loc['BST2 '+bieb[i][2], 'r2 '+o]=r_squared
                                rates.loc['BST2 '+bieb[i][2], 'Clast '+o]=name1[r].values[-1]
                                rates.loc['BST2 '+bieb[i][2], 'krel '+o]=np.absolute(sigma_ab[0]/popt[0])
                #BST4
                #Initial guess
                if name2[r].max()>0.1:
                    k0=0.001
                    A0=name2[r].values[0]
                    name2['days']=name2['hours']/24
                    popt, pcov =optimization.curve_fit(func, name2['days'], name2[r])
                    #95% confidence intervals
                    sigma_ab = np.sqrt(np.diagonal(pcov)) #https://stackoverflow.com/questions/39434402/how-to-get-confidence-intervals-from-curve-fit   
                    residuals = name2[r]- func(name2['days'], *popt)
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((name2[r]-np.mean(name2[r]))**2)
                    r_squared = 1 - (ss_res / ss_tot)
                    dt50=np.log(2)/popt[0]
                    dt50_95=np.log(2)/sigma_ab[0]
                    
                    timestep=np.arange(0,50,0.1)
                    ax[1].plot(timestep, func(timestep, *popt), '--', c='b', linewidth=0.7, label=str('%4.2f' %(popt[0])+'d$^\mathrm{-1}$'))#+', r$^\mathrm{2}$= '+('%1.2f' %(r_squared)))#+'    +-'+str('%.1E' %(sigma_ab[0]))) #'\nDT50='+str('%4.1f' %(np.log(2)/popt[0]))+'d') 
       
                    #bound_upper = func(name2['days'], (popt[0]+sigma_ab[0]), (popt[1]-sigma_ab[1]))
                    #bound_lower = func(name2['days'], (popt[0]-sigma_ab[0]), (popt[1]+sigma_ab[1]))
                    # plotting the confidence intervals
                    #ax[1].fill_between(name2['days'], bound_lower, bound_upper,
                    #                 color = 'black', alpha = 0.15)
                                   
                    ax[1].legend(loc=1,handletextpad=0.1, fontsize=5, frameon=True, facecolor='white', framealpha=1) # bbox_to_anchor=(1, 0.5),handletextpad=0.1, fontsize=5, frameon=False)
    
                    if 'Compound'in r:
                        for i in bieb:
                            if i in r:
                                rates.loc['BST4 '+bieb[i][2],'DT50 '+ o]=dt50
                                rates.loc['BST4 '+bieb[i][2],'DT50 95% '+ o]=dt50_95
                                rates.loc['BST4 '+bieb[i][2], 'k '+ o]=popt[0]
                                rates.loc['BST4 '+bieb[i][2], 'A0 '+ o]=popt[1]
                                rates.loc['BST4 '+bieb[i][2], 'k95% '+o]=sigma_ab[0]  
                                rates.loc['BST4 '+bieb[i][2], 'A095% '+o]=sigma_ab[1] 
                                rates.loc['BST4 '+bieb[i][2], 'r2 '+o]=r_squared
                                rates.loc['BST4 '+bieb[i][2], 'Clast '+o]=name2[r].values[-1]
                                rates.loc['BST4 '+bieb[i][2], 'krel '+o]=np.absolute(sigma_ab[0]/popt[0])

                #BST5
                #Initial guess                                    
                if name3[r].max()>0.1:
                    k0=0.01
                    A0=name3[r].values[0]
                    name3['days']=name3['hours']/24
                    try:
                        popt, pcov =optimization.curve_fit(func, name3['days'], name3[r])
                        #95% confidence intervals
                        sigma_ab = np.sqrt(np.diagonal(pcov)) #https://stackoverflow.com/questions/39434402/how-to-get-confidence-intervals-from-curve-fit   
                        residuals = name3[r]- func(name3['days'], *popt)
                        ss_res = np.sum(residuals**2)
                        ss_tot = np.sum((name3[r]-np.mean(name3[r]))**2)
                        r_squared = 1 - (ss_res / ss_tot)
                        dt50=np.log(2)/popt[0]
                        dt50_95=np.log(2)/sigma_ab[0]
                    
                        timestep=np.arange(0,50,0.1)
                        ax[1].plot(timestep, func(timestep, *popt), '--', c='r', linewidth=0.7, label=str('%4.2f' %(popt[0])+'d$^\mathrm{-1}$'))#+', r$^\mathrm{2}$= '+('%1.2f' %(r_squared)))#+'    +-'+str('%.1E' %(sigma_ab[0]))) #'\nDT50='+str('%4.1f' %(np.log(2)/popt[0]))+'d') 
        
                        #bound_upper = func(name3['days'], (popt[0]+sigma_ab[0]), (popt[1]-sigma_ab[1]))
                        #bound_lower = func(name3['days'], (popt[0]-sigma_ab[0]), (popt[1]+sigma_ab[1]))
                        # plotting the confidence intervals
                        #ax[1].fill_between(name3['days'], bound_lower, bound_upper,
                        #                 color = 'black', alpha = 0.15)
                                       
                        ax[1].legend(loc=1,handletextpad=0.1, fontsize=5, frameon=True, facecolor='white', framealpha=1) # bbox_to_anchor=(1, 0.5),handletextpad=0.1, fontsize=5, frameon=False)

                    except:
                        print(o+' te weinig datapunten')
                        
                    if 'Compound'in r:
                        for i in bieb:
                            if i in r:
                                rates.loc['BST5 '+bieb[i][2],'DT50 '+ o]=dt50
                                rates.loc['BST5 '+bieb[i][2],'DT50 95% '+ o]=dt50_95
                                rates.loc['BST5 '+bieb[i][2], 'k '+ o]=popt[0]
                                rates.loc['BST5 '+bieb[i][2], 'A0 '+ o]=popt[1]
                                rates.loc['BST5 '+bieb[i][2], 'k95% '+o]=sigma_ab[0]  
                                rates.loc['BST5 '+bieb[i][2], 'A095% '+o]=sigma_ab[1] 
                                rates.loc['BST5 '+bieb[i][2], 'r2 '+o]=r_squared
                                rates.loc['BST5 '+bieb[i][2], 'Clast '+o]=name2[r].values[-1]
                                rates.loc['BST5 '+bieb[i][2], 'krel '+o]=np.absolute(sigma_ab[0]/popt[0])

                ax[1].set_ylim(ymin=0)
                ax[1].grid()
                
                ax2 = ax[1].twinx()
                ax2.get_yaxis().set_visible(False)
                
                lines = ax[1].get_lines()
                if o == 'MW 5':
                    legend2 = plt.legend([lines[i] for i in [0,1,2,3,4,5]], ['Storage period - Winter 2019', 'Storage period - Fall 2020','Storage period - Spring 2021','$k_1$', '$k_2$', '$k_3$'], loc='upper center', bbox_to_anchor=(0.5, -0.4), ncol=2)
                    ax2.add_artist(legend2)

            else:
                print('Not enoug measurements for fitting '+str(r[13:-2]))
                        
        plt.subplots_adjust(hspace=0.3, wspace=0.3) 
        fig.text(0.4, -0.00, 'Time after first sample (days)') 
        
        #Add legend for all lines:        
#        ax2 = ax[1].twinx()
#        ax2.get_yaxis().set_visible(False)
        
#        lines = ax[1].get_lines()
#        if len(lines)==4:
#            legend2 = plt.legend([lines[i] for i in [0,1,2,3]], ['Storage period - Winter 2019', 'Storage period - Fall 2020','DT50', 'DT50'], loc='upper center', bbox_to_anchor=(-1.3, -0.6), ncol=2)
#            ax2.add_artist(legend2)
        
#        for ax in fig.axes:
#            ax.set_ylim(0,5)

        fig.text(0.05, 0.6, 'Concentration (mg/l)', rotation=90)
        plt.savefig('../Figures/Figure_4.png',bbox_inches='tight', dpi=300)    

#rates.to_csv('C:\\Users\\emiel.kruisdijk\\Documents\\AGRIMAR\\Work\\20170703 Mon. MAR\Breakthrough-storage test\\Overkoepelend\\Figures\\First-order\\v2 flexible A0\\20200709 Rates first-order NO3.csv')
