# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 16:29:43 2019

@author: emiel.kruisdijk
"""
import numpy as np
from matplotlib import pyplot as plt
import scipy
from pylab import *
import pandas as pd
import matplotlib as mpl
import matplotlib.patches as patches
from matplotlib.offsetbox import AnchoredText
import datetime
import functools
from ftplib import FTP
import os
import matplotlib.dates as mdates
    
def NMCpro_data():
    ########################################################################
    workspace='../Data/ASR_operation_data/'
    ########################################################################
    
    ########################################################################
    # Make Dataframe with all the data
    
    frames = []
    
    for file in os.listdir(workspace):
        if 'History Irrigation Infiltrat'in file:
            statinfo = os.stat(workspace+file)
            size=statinfo.st_size
            print(file, size)
            if size > 0:
                er=pd.read_csv(workspace+file, sep=',')
    
                frames.append(er)
            else:
                continue
    Infil = pd.concat(frames)
    Infil = Infil[['Datum', 'Water', 'Tijdsduur', 'Flow']]
    Infil['Datum']=pd.to_datetime(Infil['Datum'], format='%Y-%m-%d %H:%M:%S')
    Infil['Water']=pd.to_numeric(Infil['Water'])
    Infil['Flow']=pd.to_numeric(Infil['Flow'])
    Infil=Infil.set_index(Infil['Datum'])
    
    Infil['Water']=np.where(Infil.index < pd.to_datetime('2020-12-12 00:00:00', format='%Y-%m-%d %H:%M:%S'), Infil['Water']/1000, Infil['Water'])
#    Infil.loc[Infil.index < pd.to_datetime('2020-12-12 00:00:00', format='%Y-%m-%d %H:%M:%S'), 'Water']=Infil['Water']/1000
    Infil=Infil.sort_values(by='Datum')  
    Infil=Infil.drop_duplicates() 
    frames = []
    
    for file in os.listdir(workspace):
        if 'History Irrigation onttrek'in file:
            statinfo = os.stat(workspace+file)
            size=statinfo.st_size
            print(file, size)
            if size > 0:
                er=pd.read_csv(workspace+file, sep=',')
    
                frames.append(er)
            else:
                continue
    Onttrek = pd.concat(frames)
    Onttrek = Onttrek[['Datum', 'Water', 'Tijdsduur', 'Flow']]
    Onttrek['Datum']=pd.to_datetime(Onttrek['Datum'], format='%Y-%m-%d %H:%M:%S')
    Onttrek['Water']=pd.to_numeric(Onttrek['Water'])
    Onttrek['Flow']=pd.to_numeric(Onttrek['Flow'])
    Onttrek=Onttrek.set_index(Onttrek['Datum'])

    return Infil, Onttrek

Infil, Onttrek = NMCpro_data()
       
def CR1000_figure_period1(start_date, end_date, Infil_daily, Lt):  
    mask = (Infil_daily.index > start_date) & (Infil_daily.index <= end_date)
    Infil_daily_1=Infil_daily.loc[mask]
   
    #ax1.plot(Infil_daily_1.index, Infil_daily_1['Water_cumsum'], '.', ms=5, label='Injected water (m3)')
    ax1.plot([pd.to_datetime('2019-11-01 00:00:00', format='%Y-%m-%d %H:%M:%S'), pd.to_datetime('2019-12-10 00:00:00', format='%Y-%m-%d %H:%M:%S'), pd.to_datetime('2019-12-18 00:00:00', format='%Y-%m-%d %H:%M:%S')], [0,1374, 1374], '--')
    ax1.set_ylabel('Injected Water (m$^\mathrm{3}$)')
    ax1.set_xlim(start_date, end_date)
    ax1.spines['right'].set_visible(False)
    d = .05  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-d, -1), (d, 1)], markersize=20,
                          linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([1, 1], [0,1], transform=ax1.transAxes, **kwargs)
    ax1.xaxis.set_tick_params(which='both',bottom='on', top='off', labelbottom='on', labelrotation=45)
    ax1.set_ylim(0,5000)
    ax1.plot(Lt['Date'], Lt['Inf.1'], 'o', c='r', markersize=4, zorder=10)
    ax1.locator_params(axis='x', nbins=5)
    ax1.set_title('Operation period 1', fontsize=8.5,fontweight='bold')
    
    st=pd.to_datetime('19-12-2019', format='%d-%m-%Y')
    en=pd.to_datetime('02-02-2020', format='%d-%m-%Y')
    
    lu=[ax1]
    
    for k in lu:    
        rect = patches.Rectangle((mdates.date2num(st),-100000),(mdates.date2num(en)-mdates.date2num(st)),300000,linewidth=1,edgecolor='g',facecolor='g', alpha=0.1)
        k.add_patch(rect)   
    
    plt.locator_params(axis='x', nbins=5)
    plt.locator_params(axis='y', nbins=5)
    #f.suptitle('Infiltration water')
    
    #plt.savefig('C:/Users/emiel.kruisdijk/Documents/AGRIMAR/Work/20170703 Mon. MAR/Data/Data_check/Archive/CR1000 '+ str(start_date)[:-9] +'_'+ str(end_date)[:-9]+'.png', bbox_inches='tight', dpi=300)

def CR1000_figure_period2(start_date, end_date, Infil_daily, Lt, Infil):
   
    mask = (Infil_daily.index > start_date) & (Infil_daily.index <= end_date)
    Infil_daily_2=Infil_daily.loc[mask]

    mask = (Infil['Datum'] > start_date) & (Infil['Datum'] <= end_date)
    Nmc_2=Infil.loc[mask]    
    #Visualize data (Battery volt, temperatuur, EC, pH, eH, O2, Turb, Precipitation)
       
    ax2.plot(Infil_daily_2.index, Infil_daily_2['Water_cumsum'], '.', ms=5, label='Injected water (m$^\mathrm{3}$)')
    #ax11.set_ylabel('Injected Water\n (m3)')
    ax2.set_xlim(start_date, end_date)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    d = .05  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-d, -1), (d, 1)], markersize=20,
                          linestyle="none", color='k', mec='k', mew=1, clip_on=False)   
    ax2.plot([1, 1], [0,1], transform=ax2.transAxes, **kwargs)
    ax2.set_zorder(1)
    ax2.plot([0, 0], [1,0], transform=ax2.transAxes, **kwargs)
    ax2.set_zorder(1)    
    ax2.xaxis.set_tick_params(which='both',bottom='on', top='off', labelbottom='on', labelrotation=45)
    ax2.yaxis.set_tick_params(which='both',left='off', right='off', labelleft='off')
    ax2.set_ylim(0,5000)
    ax2.plot(Lt['Date'], Lt['Inf.1'], 'o', c='r', markersize=4, zorder=10)
    ax2.locator_params(axis='x', nbins=5)
    ax2.set_title('Operation period 2', fontsize=8.5,fontweight='bold')
    
    lu=[ax2]
    
    for k in lu:
        for index, row in Nmc_2.iterrows():
            if row['Tijdsduur'] > pd.to_timedelta('15 min'):
                rect = patches.Rectangle((mdates.date2num(row['Datum']),-100000),(mdates.date2num(row['Einddatum'])-mdates.date2num(row['Datum'])),300000,linewidth=1,edgecolor='k',facecolor='k', alpha=0.1)
                k.add_patch(rect)
        
        st=pd.to_datetime('14-10-2020', format='%d-%m-%Y')
        en=pd.to_datetime('25-11-2020', format='%d-%m-%Y')
        
        rect = patches.Rectangle((mdates.date2num(st),-100000),(mdates.date2num(en)-mdates.date2num(st)),300000,linewidth=1,edgecolor='g',facecolor='g', alpha=0.1)
        k.add_patch(rect)        

def CR1000_figure_period3(start_date, end_date, Infil_daily, Lt, Infil):   
    mask = (Infil_daily.index > start_date) & (Infil_daily.index <= end_date)
    Infil_daily_3=Infil_daily.loc[mask]
    
    mask = (Infil['Datum'] > start_date) & (Infil['Datum'] <= end_date)
    Nmc_3=Infil.loc[mask]
        
    ax3.plot(Infil_daily_3.index, Infil_daily_3['Water_cumsum'], '.', ms=5, label='Injected water (m3)')
    #ax12.set_ylabel('Injected Water\n (m3)')
    ax3.set_xlim(start_date, end_date)
    ax3.spines['left'].set_visible(False)
    d = .05  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-d, -1), (d, 1)], markersize=20,
                          linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax3.plot([0, 0], [1,0], transform=ax3.transAxes, **kwargs)
    ax3.xaxis.set_tick_params(which='both',bottom='on', top='off', labelbottom='on', labelrotation=45)
    ax3.yaxis.set_tick_params(which='both',left='off', right='off', labelleft='off')
    ax3.set_ylim(0,5000)
    ax3.plot(Lt['Date'], Lt['Inf.1'], 'o', c='r', markersize=4, zorder=10)
    ax3.locator_params(axis='x', nbins=5)
    ax3.set_title('Operation period 3', fontsize=8.5,fontweight='bold')

    lu=[ax3]
    
    for k in lu:
        for index, row in Nmc_3.iterrows():
            if row['Tijdsduur'] > pd.to_timedelta('15 min'):
                rect = patches.Rectangle((mdates.date2num(row['Datum']),-100000),(mdates.date2num(row['Einddatum'])-mdates.date2num(row['Datum'])),300000,linewidth=1,edgecolor='k',facecolor='k', alpha=0.1)
                k.add_patch(rect)
       
def Litertellers():
    ########################################################################
    workspace='../Data/'
    ########################################################################
    Lt=pd.read_csv(workspace+'Manual_readings_volumemeters.csv', sep=';', header=0)
    Lt['Date'] =  pd.to_datetime(Lt['Date'], format='%d-%m-%Y  %H:%M')
    
    return Lt
    
###############################################################################
########################## Figure  ############################################
    
########################################################################
#STYLE

#Legend
mpl.rcParams['legend.numpoints']=1
mpl.rcParams['legend.handlelength']=2
#mpl.legend(loc='best', frameon=True)
mpl.rcParams['legend.frameon'] = True
mpl.rcParams['legend.loc'] = 'best'
mpl.rcParams['legend.fontsize'] = 11
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.markersize'] = 1
#######################################################################
#Visualize data (Battery volt, temperatuur, EC, pH, eH, O2, Turb, Precipitation)

f, ((ax1, ax2, ax3))= subplots(nrows=1, ncols=3, figsize=(6,1.8), gridspec_kw={'width_ratios': [1,1,1]})

###############################################################################
################### NMCpro data updaten #######################################

Infil, Onttrek = NMCpro_data()

###############################################################################
##################### Add data from first period infiltration #################
###############################################################################
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

df1=System[['System_start', 'Volume_start_INF1', 'Volume_start_INF2']]
df1.columns= ['Datum', 'INF1', 'INF2']
df2=System[['System_end', 'Volume_end_INF1', 'Volume_end_INF2']]
df2.columns= ['Datum', 'INF1', 'INF2']
df3=pd.concat([df1, df2])
df3=df3.sort_values(by=['Datum'])
df3['Water'] = df3['INF1'].shift(1).fillna(0)
df3['Water'] = df3['INF1'] - df3['Water']
df3=df3.drop(['INF1','INF2'], axis=1)
df3['Tijdsduur']=0
df3['Flow']=0

#add row just before stop infiltration:
listd= [[pd.to_datetime('2019-12-04 10:08:00', format='%Y-%m-%d %H:%M:%S'), 431.06, 0, 0], [pd.to_datetime('2019-12-10 00:00:00', format='%Y-%m-%d %H:%M:%S'), 471.13, 0, 0]]  #Water is 902.46 m3 gemeten op 4 december - totaal tijdens BST1 471.4 = 431.06

df3=df3.append(pd.DataFrame(listd, columns=['Datum', 'Water', 'Tijdsduur', 'Flow']))
df3=df3.set_index(df3['Datum'])

Infil=pd.concat([df3, Infil])

Infil = Infil[(Infil.index != pd.to_datetime('30-10-2020  00:00:00', format='%d-%m-%Y %H:%M:%S')) & (Infil.index != pd.to_datetime('29-10-2020  08:35:00', format='%d-%m-%Y %H:%M:%S'))]

Infil['Tijdsduur']=pd.to_timedelta(Infil['Tijdsduur'])

#Infil = Infil[(Infil['Tijdsduur']== pd.to_timedelta('00:00:00')) | (Infil['Tijdsduur']>pd.to_timedelta('00:01:00'))]

Infil = Infil.loc[Infil['Water']!=0]
Infil['Einddatum']=Infil['Datum']+Infil['Tijdsduur']

Infil_daily=Infil.resample('D').agg({'Water': np.sum, 'Flow': np.mean})

Infil_daily['Water_cumsum']= Infil_daily['Water'].cumsum()


########################## PERIOD 1 ###########################################
    
start_date=pd.to_datetime('01-10-2019', format='%d-%m-%Y')
end_date=pd.to_datetime('19-12-2019', format='%d-%m-%Y')

Lt=Litertellers()

CR1000_figure_period1(start_date, end_date, Infil_daily,  Lt)

###############################################################################
########################## PERIOD 2 ###########################################

start_date=pd.to_datetime('19-09-2020', format='%d-%m-%Y')
end_date=pd.to_datetime('14-10-2020', format='%d-%m-%Y')

CR1000_figure_period2(start_date, end_date, Infil_daily,  Lt, Infil)

###############################################################################
########################## PERIOD 3 ###########################################

start_date=pd.to_datetime('05-12-2020', format='%d-%m-%Y')
end_date= pd.to_datetime('18-03-2021', format='%d-%m-%Y')  

CR1000_figure_period3(start_date, end_date, Infil_daily, Lt, Infil)

###############################################################################
################# ADD INFORMATIONAL TEXT TO FIGURE ############################
f.text(0.37, 0.75, 'Storage period 1', fontsize=9, rotation='vertical')
f.text(0.65, 0.75, 'Storage period 2', fontsize=9, rotation='vertical')
f.text(0.915, 0.75, 'Storage period 3', fontsize=9, rotation='vertical')

ax1.arrow(pd.to_datetime('21-10-2019', format='%d-%m-%Y'), 3500, -12, 0, head_width=400, head_length=5, edgecolor='k', facecolor='k')
ax1.text(pd.to_datetime('10-10-2019', format='%d-%m-%Y'), 3700, 'Push-Pull \n Tests')
ax1.arrow(pd.to_datetime('31-10-2019', format='%d-%m-%Y'), 1500, 0, -1300, head_width=7.5, head_length=200, edgecolor='k', facecolor='k')
ax1.text(pd.to_datetime('07-10-2019', format='%d-%m-%Y'), 2000, 'Continuous Field \nInjection Experiment')

plt.subplots_adjust(hspace=0.2, wspace=0.3) 
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)
plt.savefig('../Figures/Figure_1.png',bbox_inches='tight', dpi=300)    
