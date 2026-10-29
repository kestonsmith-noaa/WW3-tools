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

#import itertools

#import FiniteElementMeshRoutines as FE
#import GaussMarkov as GM

print("Defining meshname and other variables from file DefineVariables.py")
with open("DefineVariables.py", "r") as file:
    code = file.read()
    print(code)
    exec(code)

xs=[]
ys=[]
zs=[]
es=[]
#es=np.array(es)
#for j in range(nf):
j=nf-1
fl=flnms[j]
data = nc.Dataset(BathyDir+fl,"r")
x0=data["lon"][:]
y0=data["lat"][:]
print(fl)
skip=3
skip=10 #Tiny
#skip=20 #Nano
if filetype[j] <2 : #gridded data
    Z0=data["z"][:,:]
    nx=len(x0)
    ny=len(y0)
    X0=np.ones((ny,1))*np.transpose(x0[:])
    Y0=np.transpose(np.ones((nx,1))*np.transpose(y0[:]))
    X0=X0[0 :: skip,0 :: skip]
    Y0=Y0[0 :: skip,0 :: skip]
    Z0=Z0[0 :: skip,0 :: skip]
    x0 = X0.flatten() 
    y0 = Y0.flatten() 
    z0 = Z0.flatten() 
else :
    z0=data["z"][:]
        
fillval0=data["z"]._FillValue
jf=np.where( z0 != fillval0 )
x0=np.array(x0[jf])
y0=np.array(y0[jf])
z0=np.array(z0[jf])

jf=np.where( z0**2 >= 0. )
x0=np.array(x0[jf])
y0=np.array(y0[jf])
z0=np.array(z0[jf])
    
es=np.concatenate((es, z0*0. + fileStdErr[j]))
xs=np.concatenate((xs, x0))
ys=np.concatenate((ys, y0))
zs=np.concatenate((zs, z0))

TextOutput=False
if TextOutput :
    np.savetxt("GlobalDataBathyX.txt" , xs, fmt="%.8e")
    np.savetxt("GlobalDataBathyY.txt" , ys, fmt="%.8e" )
    np.savetxt("GlobalDataBathyZ.txt" , zs, fmt="%.8e" )
    np.savetxt("GlobalDataBathyStd.txt" , es, fmt="%.8e" )
else :    
#   with nc.Dataset("GlobalDataPointsTinySmall.nc", 'w', format='NETCDF4') as ncout:
   with nc.Dataset("GlobalDataPointsTiny.nc", 'w', format='NETCDF4') as ncout:
        n=len(zs)
        ncout.createDimension('N', n)  # Unlimited dimension
        
        lon_var=ncout.createVariable('lon', 'f4', ('N',))
        lon_var.units         = 'degree_east'
        lon_var.long_name     = 'longitude'
        lon_var.standard_name = 'longitude'
        lon_var.axis          = 'lon'
        lon_var[:]=xs[:]

        lat_var=ncout.createVariable('lat', 'f4', ('N',))
        lat_var.units         = 'degree_north'
        lat_var.long_name     = 'latitude'
        lat_var.standard_name = 'latitude'
        lat_var.axis          = 'lat'
        lat_var[:]=ys[:]

        z_var=ncout.createVariable('z', 'f4', ('N',))
        z_var.long_name     = 'z'
        z_var.units         = 'meters'
        z_var.positive      = 'up'
        z_var.standard_name = 'height'
        z_var[:]=zs[:]

        e_var=ncout.createVariable('StdErr', 'f4', ('N',))
        e_var.long_name     = 'z'
        e_var.units         = 'meters'
        e_var.positive      = 'up'
        e_var.standard_name = 'height'
        e_var[:]=es[:]
        
        ncout.close
