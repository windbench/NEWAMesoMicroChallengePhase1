#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 09:22:33 2017

@author: usuario
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import pandas as pd
import netCDF4 
from scipy.interpolate import interp1d
import seaborn as sns
from windrose import WindroseAxes
from windrose import plot_windrose
import matplotlib.cm as cm

# Constants
g = 9.81    # [m s-2]
P0 = 100000 # Reference pressure [Pa]
T0 = 300    # Reference temperature for perturbation temperature [K]
kappa = 0.2854  # Poisson constant (R/Cp)
R_air = 287.058  # Specific gas constant for dry air [J kg-1 K-1]
Cp_air = 1005   # Specific heat of air [J kg-1 K-1]
omega = 7.2921159e-5    # angular speed of the Earth [rad/s]
K = 0.41    # von Karman constant

# Site 
siteID = 'Cabauw'
lat_s = 51.971   # degrees N
lon_s = 4.927    # degrees E
fc  = 2*omega*np.sin(lat_s*np.pi/180)     # Coriolis parameter [s-1]
datadescription = 'Profile data from the CESAR observatory at KNMI (http://www.cesar-observatory.nl/)'

# Evaluation period
datefrom = datetime.datetime(2006,1,1,0,0,0)
dateto = datetime.datetime(2006,12,30,23,50,0)
savenc = 1

# Load data from CESAR
dirobs = './CESAR/'

filesufix = []
for ye in range(datefrom.year, dateto.year+1):
    if ye == datefrom.year:
        for mo in range(datefrom.month, 13):
            filesufix.append('{0:d}{1:02d}'.format(ye,mo))
    elif ye == dateto.year:
        for mo in range(1, dateto.month +1):
            filesufix.append('{0:d}{1:02d}'.format(ye,mo))
    else:
        for mo in range(1, 13):
            filesufix.append('{0:d}{1:02d}'.format(ye,mo))
            
nodata = -9999.0    # missing data flag

tmet_obs = [];
WD10_obs = []; S10_obs = []; P0_obs = []; RH2_obs = []; T2_obs = [];

# Surface met data
fileobs = 'cesar_surface_meteo_lc1_t10_v1.0_'
for f in range (0,len(filesufix)):
    file1 = netCDF4.Dataset(dirobs+fileobs+filesufix[f]+'.nc', 'r')
    
    #f.variables['ustar'].long_name   # name of variable
    
    date_obs = file1.variables['date'][:]
    time_obs = file1.variables['time'][:]   # 'hours since 2015-01-01 00:00:00 (UTC)'   
    time_obs = time_obs * 3600          
    datetime_obs = []
    date0 = datetime.datetime(int(filesufix[f][0:4]),int(filesufix[f][4:6]),1,0,0,0)                   
    for i in range (0,len(time_obs)):
        datetime_obs.append(date0 + datetime.timedelta(seconds = np.int(time_obs[i])))
    
    ifrom_obs=0
    for j in range(1,len(datetime_obs)):
        if datetime_obs[j] <= datefrom:
            ifrom_obs = j
    ito_obs=0
    for j in range(1,len(datetime_obs)):
        if datetime_obs[j] <= dateto:
            ito_obs = j
    
    tmet_obs = np.append(tmet_obs, mdates.date2num(datetime_obs[ifrom_obs:ito_obs]))
     
    # Time series
    WD10_obs = np.append(WD10_obs, np.ma.array(file1.variables['D010'][ifrom_obs:ito_obs], mask=(file1.variables['D010'][ifrom_obs:ito_obs] == nodata)))
    S10_obs  = np.append(S10_obs, np.ma.array(file1.variables['F010'][ifrom_obs:ito_obs], mask=(file1.variables['F010'][ifrom_obs:ito_obs] == nodata)))
    P0_obs   = np.append(P0_obs, np.ma.array(file1.variables['P0'][ifrom_obs:ito_obs], mask=(file1.variables['P0'][ifrom_obs:ito_obs] == nodata)))
    RH2_obs  = np.append(RH2_obs, np.ma.array(file1.variables['RH002'][ifrom_obs:ito_obs], mask=(file1.variables['RH002'][ifrom_obs:ito_obs] == nodata)))   # relative humidity
    T2_obs   = np.append(T2_obs, np.ma.array(file1.variables['TA002'][ifrom_obs:ito_obs], mask=(file1.variables['TA002'][ifrom_obs:ito_obs] == nodata)))    # relative humidity

T2_obs = T2_obs + 273.15

# Surface flux data
fileobs = 'cesar_surface_flux_lc1_t10_v1.0_'
tflux_obs = []
us_obs = []; H_obs = []; G0_obs = []
for f in range (0,len(filesufix)):
    file1 = netCDF4.Dataset(dirobs+fileobs+filesufix[f]+'.nc', 'r')
    
    #file1.variables['ustar'].long_name   # name of variable
    
    date_obs = file1.variables['date'][:]
    time_obs = file1.variables['time'][:]   # 'hours since 2015-01-01 00:00:00 (UTC)'   
    time_obs = time_obs * 3600          
    datetime_obs = []
    date0 = datetime.datetime(int(filesufix[f][0:4]),int(filesufix[f][4:6]),1,0,0,0)                   
    for i in range (0,len(time_obs)):
        datetime_obs.append(date0 + datetime.timedelta(seconds = np.int(time_obs[i])))
    
    ifrom_obs=0
    for j in range(1,len(datetime_obs)):
        if datetime_obs[j] <= datefrom:
            ifrom_obs = j
    ito_obs=0
    for j in range(1,len(datetime_obs)):
        if datetime_obs[j] <= dateto:
            ito_obs = j
    
    tflux_obs = np.append(tflux_obs, mdates.date2num(datetime_obs[ifrom_obs:ito_obs]))
     
    # Time series
    us_obs = np.append(us_obs, np.ma.array(file1.variables['UST'][ifrom_obs:ito_obs], mask=(file1.variables['UST'][ifrom_obs:ito_obs] == nodata)))
    H_obs =  np.append(H_obs, np.ma.array(file1.variables['H'][ifrom_obs:ito_obs], mask=(file1.variables['H'][ifrom_obs:ito_obs] == nodata)))            # Surface sensible heat flux 
    G0_obs = np.append(G0_obs, np.ma.array(file1.variables['G0'][ifrom_obs:ito_obs], mask=(file1.variables['G0'][ifrom_obs:ito_obs] == nodata)))      # Surface soil heat flux

zflux = 3.

rho_obs = P0_obs*100./(R_air*(T2_obs))
wt_obs = H_obs/(rho_obs*Cp_air)


# Tower met data
fileobs = 'cesar_tower_meteo_lb1_t10_v1.2_'
tmet_obs = []
WD_obs = []; S_obs = []; Sstd_obs = []; T_obs = []; RH_obs = []; 
for f in range (0,len(filesufix)):
    file1 = netCDF4.Dataset(dirobs+fileobs+filesufix[f]+'.nc', 'r')

    date_obs = file1.variables['date'][:]
    time_obs = file1.variables['time'][:]   # 'hours since 2015-01-01 00:00:00 (UTC)'   
    time_obs = time_obs * 3600          
    datetime_obs = []
    date0 = datetime.datetime(int(filesufix[f][0:4]),int(filesufix[f][4:6]),1,0,0,0)                   
    for i in range (0,len(time_obs)):
        datetime_obs.append(date0 + datetime.timedelta(seconds = np.int(time_obs[i])))
    
    ifrom_obs=0
    for j in range(1,len(datetime_obs)):
        if datetime_obs[j] <= datefrom:
            ifrom_obs = j
    ito_obs=0
    for j in range(1,len(datetime_obs)):
        if datetime_obs[j] <= dateto:
            ito_obs = j
    
    tmet_obs = np.append(tmet_obs, mdates.date2num(datetime_obs[ifrom_obs:ito_obs]))
         
    # Time series
    z_obs = file1.variables['z'][ifrom_obs:ito_obs]   # tower heights 
    WD_obs = np.ma.append(WD_obs, file1.variables['D'][ifrom_obs:ito_obs])
    S_obs = np.ma.append(S_obs, file1.variables['F'][ifrom_obs:ito_obs])
    Sstd_obs = np.ma.append(Sstd_obs, file1.variables['SF'][ifrom_obs:ito_obs])
#    RH_obs = np.ma.append(RH_obs, file1.variables['RH'][ifrom_obs:ito_obs])
    T_obs = np.ma.append(T_obs, file1.variables['TA'][ifrom_obs:ito_obs])

S_obs = np.reshape(S_obs,(-1,len(z_obs)))
Sstd_obs = np.reshape(Sstd_obs,(-1,len(z_obs)))
WD_obs = np.reshape(WD_obs,(-1,len(z_obs)))
#RH_obs = np.reshape(RH_obs,(-1,len(z_obs)))
T_obs = np.reshape(T_obs,(-1,len(z_obs))) 

T0_obs = T2_obs 
rho0_obs = P0_obs*1e2/(R_air*T0_obs)
wt0_obs = H_obs/(rho0_obs*Cp_air)
L_obs = -us_obs**zflux/(K*(g/T0_obs)*wt0_obs)
L_obs[np.where(abs(L_obs)<0.1)] = np.nan     # Remove outliers when L is too small
zL_obs = zflux/L_obs

WD10_obs = pd.DataFrame(WD10_obs, index = mdates.num2date(tmet_obs))
S10_obs = pd.DataFrame(S10_obs, index = mdates.num2date(tmet_obs))
P0_obs = pd.DataFrame(P0_obs, index = mdates.num2date(tmet_obs))
RH2_obs = pd.DataFrame(RH2_obs, index = mdates.num2date(tmet_obs))
T2_obs = pd.DataFrame(T2_obs, index = mdates.num2date(tmet_obs))

us_obs = pd.DataFrame(us_obs, index = mdates.num2date(tflux_obs))
H_obs = pd.DataFrame(H_obs, index = mdates.num2date(tflux_obs))
G0_obs = pd.DataFrame(G0_obs, index = mdates.num2date(tflux_obs))
wt_obs = pd.DataFrame(wt_obs, index = mdates.num2date(tflux_obs))
L_obs = pd.DataFrame(L_obs, index = mdates.num2date(tflux_obs))

S_obs = pd.DataFrame(S_obs, index = mdates.num2date(tmet_obs), columns = z_obs)
Sstd_obs = pd.DataFrame(Sstd_obs, index = mdates.num2date(tmet_obs), columns = z_obs)
WD_obs = pd.DataFrame(WD_obs, index = mdates.num2date(tmet_obs), columns = z_obs)
#RH_obs = pd.DataFrame(RH_obs, index = mdates.num2date(tmet_obs), columns = z_obs)
T_obs = pd.DataFrame(T_obs, index = mdates.num2date(tmet_obs), columns = z_obs) 

Nt = len(S_obs)

ts = (S_obs.index.values - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
dates = [datetime.datetime.utcfromtimestamp(x) for x in ts]

# save to netcdf
if savenc == 1:
    fileout = siteID + '_mast_'+ datefrom.strftime("%Y%m%d") + '_' + dateto.strftime("%Y%m%d") + '.nc'
    
    f = netCDF4.Dataset(fileout, 'w', format='NETCDF4')
    f.history = datadescription
    
    f.createDimension('time', Nt)
    f.createDimension('zS', )
    f.createDimension('zWD', )
    f.createDimension('zT', )
    f.createDimension('zTw', )
    f.createDimension('zflux', )
    f.createDimension('site', 1)
    
    lats = f.createVariable('lat', 'float', ('site',))
    lats.long_name = 'Site latitude'
    lats.units = 'degrees North'
    lats[:] = lat_s
    
    lons = f.createVariable('lon', 'float', ('site',))
    lons.long_name = 'Site longitude'
    lons.units = 'degrees East'
    lons[:] = lon_s
    
    fcs = f.createVariable('fc', 'float', ('site',))
    fcs.long_name = 'Coriolis parameter'
    fcs.units = 's-1'
    fcs[:] = fc
    
    timess = f.createVariable('time', 'float', ('time',))
    timess.long_name = 'Time'
    timess.units = 'hours since 0001-01-01 00:00:00.0'
    timess.calendar = 'gregorian'
    timess[:] = netCDF4.date2num(dates,units=timess.units,calendar=timess.calendar)
    
    heightss_zS = f.createVariable('zS', 'float', ('zS',))
    heightss_zS.long_name = 'Height above sea level cup anemometer'
    heightss_zS.units = 'm'
    heightss_zS[:] = z_obs
    
    heightss_zWD = f.createVariable('zWD', 'float', ('zWD',))
    heightss_zWD.long_name = 'Height above sea level wind vane'
    heightss_zWD.units = 'm'
    heightss_zWD[:] = z_obs

    heightss_zT = f.createVariable('zT', 'float', ('zT',))
    heightss_zT.long_name = 'Height above sea level air temperature'
    heightss_zT.units = 'm'
    heightss_zT[:] = z_obs
        
    heightss_zflux = f.createVariable('zflux', 'float', ('zflux',))
    heightss_zflux.long_name = 'Height above sea level sonic anemometer'
    heightss_zflux.units = 'm'
    heightss_zflux[:] = zflux
    
    Ss = f.createVariable('S', 'float', ('time','zS',))
    Ss.long_name = 'Horizontal mean wind speed' 
    Ss.units = 'm s-1'
    Ss[:] = S_obs.as_matrix()
        
    Sstds = f.createVariable('Sstd', 'float', ('time','zS',))
    Sstds.long_name = 'Standard deviation of horizontal mean wind speed' 
    Sstds.units = 'm s-1'
    Sstds[:] = Sstd_obs.as_matrix()

    WDs = f.createVariable('WD', 'float', ('time','zWD',))
    WDs.long_name = 'Wind direction' 
    WDs.units = 'clockwise degrees from North'
    WDs[:] = WD_obs.as_matrix()
            
    Ts = f.createVariable('T', 'float', ('time','zT',))
    Ts.long_name = 'Air temperature' 
    Ts.units = 'K'
    Ts[:] = T_obs.as_matrix()
    
    uss = f.createVariable('us', 'float', ('time',))
    uss.long_name = 'Friction velocity' 
    uss.units = 'm s-1'
    uss[:] = us_obs.as_matrix()
               
    wts = f.createVariable('wt', 'float', ('time',))
    wts.long_name = 'Kinematic Upward sensible heat flux at surface wt = HFX/(rho*Cp)' 
    wts.units = 'K m s-1'
    wts[:] = wt_obs.as_matrix()
    
    Ls = f.createVariable('L', 'float', ('time',))
    Ls.long_name = 'M-O length' 
    Ls.units = 'm'
    Ls[:] = L_obs.as_matrix()

    f.close()
    print 'Saving ' + fileout
    
    
    

