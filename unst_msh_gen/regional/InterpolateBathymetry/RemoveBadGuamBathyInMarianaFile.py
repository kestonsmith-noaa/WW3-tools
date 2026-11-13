######################################################################
# This is the main interpolation engine.  This script reads in a list 
# of nodes created by DivideMeshNodes.py and interpolates bathymetry 
# from a list of bathymetry data files to the nodes in the assigned
# list.  The intent is to run this in NPart parallel jobs for a 
# partition of mesh nodes into NPart sets.
#
# The user is responsible for specifing the bathymetry data sets
# in list "flnms" as well as the type of data each of these files 
# in list "filetype". filetype[k] = 0 gridded local data, 1 gridded 
# global data and 2 for scattered data. 
#
# The user should define there observational error model and bacground
# process statistical model. 
#
import os
import argparse
import time
import numpy as np
import netCDF4 as nc

import sys
from scipy.interpolate import RegularGridInterpolator

import itertools

import FiniteElementMeshRoutines as FE
import GaussMarkov as GM

print("Defining meshname and other variables from file DefineVariables.py")
with open("DefineVariables.py", "r") as file:
    code = file.read()
    print(code)
    exec(code)

fl_mariana=BathyDir+"PIX/ngdc_bathy_180m_mariana.CRM.nc"
fl_guam=BathyDir+"pibhmc_bathy_60m_guam.crm.nc.S250m.VB.nc"
flout=BathyDir+"PIX/ngdc_bathy_180m_mariana.CRM.no_guam.nc"

data0=nc.Dataset(fl_mariana,"r")
data1=nc.Dataset(fl_guam,"r")
x0=data0["lon"][:]
y0=data0["lat"][:]
x1=data1["lon"][:]
y1=data1["lat"][:]
# bad data in mariana around guam.
jE=np.where( x0 > np.min(x1) )[0].tolist()
jW=np.where( x0 < np.max(x1) )[0].tolist()
jN=np.where( y0 > np.min(y1) )[0].tolist()
jS=np.where( y0 < np.max(y1) )[0].tolist()
jxb=list( set(jE) & set(jW) )
jyb=list( set(jN) & set(jS) )
z0=data0["z"][:,:]
print(z0.shape)
print(x0.shape)
print(y0.shape)
#z0[jxb,jyb]=float("nan")
for k in range(len(jxb)):
    z0[jyb,jxb[k]]=float("nan")
nx=len(x0)
ny=len(y0)
filval=-99999
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
    ncout.createDimension('lon', nx)  # Unlimited dimension
    ncout.createDimension('lat', ny)
        
    lon_var=ncout.createVariable('lon', 'f8', ('lon',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'lon'
    lon_var[:]=x0[:]

    lat_var=ncout.createVariable('lat', 'f8', ('lat',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'lat'
    lat_var[:]=y0[:]

    z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = filval)
    z_var.grid_mapping  = 'crs'
    z_var.long_name     = 'z'
    z_var.units         = 'meters'
    z_var.positive      = 'up'
    z_var.standard_name = 'height'
    z_var.vert_crs_name = 'EGM2008'
    z_var. vert_crs_epsg = 'EPSG:3855'
    z_var[:,:]=z0[:,:]
        
    ncout.close

