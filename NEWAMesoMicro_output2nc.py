# -*- coding: utf-8 -*-
"""
@author: Javier Sanz Rodrigo, jsrodrigo@cener.com, CENER
1 March 2018
"""
import numpy as np
import netCDF4 
import matplotlib.dates as mdates
import datetime

# Constants
P0      = 100000          # Reference pressure [Pa]
g       = 9.81            # [m s-2]
kappa   = 0.2854          # Poisson constant (R/Cp)
R_air   = 287.058         # Specific gas constant for dry air [J kg-1 K-1]
Cp_air  = 1005            # Specific heat of air [J kg-1 K-1]
omega   = 7.2921159e-5    # angular speed of the Earth [rad/s]
K       = 0.41            # von Karman constant

siteID = 'Cabauw'
lat_s = 51.971   # degrees N
lon_s = 4.927    # degrees E
datefrom = datetime.datetime(2006,1,1,0,0,0)
dateto = datetime.datetime(2006,12,30,23,50,0)

# Simulation description
# Please provide contact information and simulation ID and description.  
usrID = 'Javier Sanz Rodrigo (CENER)' # contact person and affiliation
emailID = 'jsrodrigo@cener.com'       # e-mail of contact person
simID = 'WRF-YSU'            # (short) name of the model. This will be used to identify your simulations in the benchmark results/plots 
simdes = 'WRFv3.8, ERA Interim, YSU, 27>9>3 km, profiles at Cabauw for 2006 from D03 spatially-averaged over L=9km'         
        # (long) description of the simulation with important config settings

# Load your simulation data 
dirsim = '.'
filesim = dirsim + '/Cabauw_tendencies_w60_L9000.nc'
f = netCDF4.Dataset(filesim, 'r')
times = f.variables['time'][:] # Days since 001-01-01 00:00:00 UTC, plus one
idates = np.where(np.logical_and(times >= mdates.date2num(datefrom), 
                             times < mdates.date2num(dateto)))[0] 
t = f.variables['time'][idates]
z = f.variables['z'][:]
U = f.variables['U'][idates,:]
V = f.variables['V'][idates,:]
Th = f.variables['Th'][idates,:]
us = f.variables['ust'][idates]
wt = f.variables['wt'][idates]
T2 = f.variables['T2'][idates]
zflux = 10.

# output variables
t_out = t
z_out = z
zflux_out = zflux
n = len(t_out)
m = len(z_out)

# time-height fields
# 2D arrays with time in rows and heights in columns using the same order of t_out and z_out
U_out = U                           # WE velocity component in natural coordinates [m s-1]
V_out = V                           # SN velocity component in natural coordinates [m s-1]
Th_out = Th                         # potential temperature
us_out = us                         # friction velocity at zflux above ground level [m s-1]
wt_out = wt                         # kinematic heat flux at zflux above ground level [K m s-1]
T2_out  = T2                        # 2-m temperature [K]
L_out = np.NAN*np.ones((n))         # Obukhov length [m]
tke_out = np.NAN*np.ones((n,m))     # turbulent kinetic energy [m2 s-2]

# Save NetCDF output file
fileout = simID + '.nc'
f = netCDF4.Dataset(fileout, 'w')
f.history = usrID+', '+emailID+', '+siteID+', '+simID+', '+simdes
f.createDimension('time', n)
f.createDimension('z', m)
f.createDimension('site', 1)

lats = f.createVariable('lat', 'float', ('site',))
lats.long_name = 'Site latitude'
lats.units = 'degrees North'
lats[:] = lat_s

lons = f.createVariable('lon', 'float', ('site',))
lons.long_name = 'Site longitude'
lons.units = 'degrees East'
lons[:] = lon_s

zfluxs = f.createVariable('zflux', 'float', ('site',))
zfluxs.long_name = 'Surface-layer flux height'
zfluxs.units = 'm'
zfluxs[:] = zflux_out
        
times = f.createVariable('time', 'float', ('time',))
times.long_name = 'Time'
times.units = 'Days since 001-01-01 00:00:00 UTC, plus one'
times[:] = t_out

heights = f.createVariable('z', 'float', ('z',))
heights.long_name = 'Height above ground level'
heights.units = 'm'
heights[:] = z_out

Us = f.createVariable('U', 'float', ('time','z',))
Us.long_name = 'U velocity component' 
Us.units = 'm s-1'
Us[:] = U_out

Vs = f.createVariable('V', 'float', ('time','z',))
Vs.long_name = 'V velocity component' 
Vs.units = 'm s-1'
Vs[:] = V_out
    
Ths = f.createVariable('Th', 'float', ('time','z',))
Ths.long_name = 'Potential temperature' 
Ths.units = 'K'
Ths[:] = Th_out

TKEs = f.createVariable('TKE', 'float', ('time','z',))
TKEs.long_name = 'Turbulent kinetic energy' 
TKEs.units = 'm2 s-2'
TKEs[:] = tke_out

uss = f.createVariable('ust', 'float', ('time',))
uss.long_name = 'Surface-layer friction velocity' 
uss.units = 'm s-1'
uss[:] = us_out

wts = f.createVariable('wt', 'float', ('time',))
wts.long_name = 'surface-layer kinematic heat flux' 
wts.units = 'K m s-1'
wts[:] = wt_out

Ls = f.createVariable('L', 'float', ('time',))
Ls.long_name = 'Obukhov length' 
Ls.units = 'm'
Ls[:] = L_out

T2s = f.createVariable('T2', 'float', ('time',))
T2s.long_name = '2-m temperature' 
T2s.units = 'K'
T2s[:] = T2_out