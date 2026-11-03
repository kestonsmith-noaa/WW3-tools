import numpy as np
import netCDF4 as nc
import scipy.sparse as sp
import os
import sys

####################################################################################################
# Converts an ascii file containing interpolation weight information for linearly interpolating from
# the nodes of an unstructured mesh to the nodes of a second different unstructed mesh into netcdf
# format. 
# Call as:
# $python convert_weights_to_netcdf.py data_file_to_be_interpolated.nc unstructured_mesh_geometry_file.msh interpolation_weights.txt extrapolation_flag
# where extrapolation_flag=True to use nearest node extrapolation, for target nodes outside the 
# coverage of the source mesh. Will create output file interpolation_weights.nc including weights 
# as well as source and destination mesh node geometry information.
#
# Input file (interpolation_weights.txt) has line format:
# column 1      : node number of target mesh
# column 2      : Source mesh element number (not used)
# columns 3,4,5 : node numbers of source mesh
# column 6      : Distance of target node to source element center (not used)
# columns 7,8,9 : Interpolation weights for interpolating from source mesh to destination 
#               : corresponding to the nodes in columns 3,4,5
####################################################################################################

nargin = len(sys.argv) - 1

flin=sys.argv[1]
mshfl=sys.argv[2]

flinslash=flin.rfind('/')+1
flinname=flin[flinslash:-3]
meshslash=mshfl.rfind('/')+1
meshname=mshfl[meshslash:-3]
flinname=flin[flinslash:-3]

TmpOutDir="InterpWeights."+meshname+flinname
print(TmpOutDir)
Njobs=int(sys.argv[3])
AllData=arr = np.zeros((1, 12))
for k in range(Njobs):
    fl=TmpOutDir+"/Part.IntrpWghts."+str(k)+".txt"
    AllOutMatrix=np.loadtxt(fl)
    print(AllOutMatrix.shape)
    FoundElement=AllOutMatrix[:,10].astype(int)
    jg=np.where(FoundElement==1)
    jg=jg[0].tolist()
    AllData=np.vstack((AllData,AllOutMatrix[jg,:]))

AllData=AllData[1 : ,:]

flout="InterpWeights."+meshname+flinname+".nc"

N,twelve=AllData.shape

with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
    ncout.createDimension('n' , N)
    ncout.createDimension('three' , 3)
    
    w_var=ncout.createVariable('weights', 'f8', ('n','three'))
    w_var.long_name     = 'node to observation interpolation weights'
    w_var[:,:]=AllData[:,6 : 9]

    n_var=ncout.createVariable('nodes', 'i4', ('n','three'))
    n_var.long_name     = 'mesh nodes'
    n_var[:,:]=AllData[:,3 : 6]
    
    en_var=ncout.createVariable('element', 'i4', ('n',))
    en_var.long_name     = 'mesh element number'
    en_var[:]=AllData[:,10]

    z_var=ncout.createVariable('z', 'f8', ('n',))
    z_var.long_name     = 'bathymetry z observed'
    z_var.units         = 'm'
    z_var[:]=AllData[:,1]

    e_var=ncout.createVariable('err', 'f4', ('n',))
    e_var.long_name     = 'standard deviation of observation error'
    e_var.units         = 'm'
    e_var[:]=AllData[:,2]
    
    idp_var=ncout.createVariable('data_point', 'i4', ('n',))
    idp_var.long_name     = 'initial data point number in :'+flin
    idp_var[:]=AllData[:,0]
    
    ncout.close


# construct sparse interpolation matrix for csr format
n_s=3*N
weights=AllData[:,6 : 9]
nodes  =AllData[:,3 : 6]
row = range(N)
ROW = np.column_stack([row] * 3) 

flout="SparseMatrix."+meshname+flinname+".nc"

nn_dst=N
import FiniteElementMeshRoutines as FE
x,y,z,e,bnd=FE.loadWW3Mesh(mshfl)
nn_src=len(x)
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
    ncout.createDimension('n_s' , n_s)
    ncout.createDimension('n' , N)
    ncout.setncattr("Nrows", nn_dst)
    ncout.setncattr("Ncols", nn_src)
    ncout.setncattr("SrcFieldType", "unstructured") 
    
    r_var=ncout.createVariable('row', 'i4', ('n_s',))
    r_var.long_name     = 'row index'
    r_var[:]=ROW.flatten()
    
    c_var=ncout.createVariable('col', 'i4', ('n_s',))
    c_var.long_name     = 'column index'
    c_var[:]=nodes.flatten()
    
    s_var=ncout.createVariable('S', 'f8', ('n_s',))
    s_var.long_name     = 'matrix value'
    s_var[:]=weights.flatten()
    
    z_var=ncout.createVariable('z', 'f8', ('n',))
    z_var.long_name     = 'bathymetry z observed'
    z_var.units         = 'm'
    z_var[:]=AllData[:,1]

    e_var=ncout.createVariable('err', 'f4', ('n',))
    e_var.long_name     = 'standard deviation of observation error'
    e_var.units         = 'm'
    e_var[:]=AllData[:,2]

    ncout.close

# Search across parallel jobs to establish local observed min, max, and mean
nn=len(x)

Zmin=np.zeros(nn)+float('inf')
Zmax=np.zeros(nn)-float('inf')
Zsum=np.zeros(nn)
Nobs=np.zeros(nn,dtype=int)

for k in range(Njobs):
    fl=TmpOutDir+"/Part.zlim."+str(k)+".txt"
    print(fl)
    AllOutMatrix=np.loadtxt(fl)
    print(AllOutMatrix.shape)

    nj0=AllOutMatrix[:,0].astype(int)
    Zmin0=AllOutMatrix[:,1]
    Zmax0=AllOutMatrix[:,2]
    Nobs0=AllOutMatrix[:,3]
    Zsum0=AllOutMatrix[:,4]

    Zmin[nj0]=np.min(np.vstack( (Zmin[nj0],Zmin0) ), axis=0 )
    Zmax[nj0]=np.max(np.vstack( (Zmax[nj0],Zmax0) ), axis=0 )
    Zsum[nj0]=Zsum[nj0]+Zsum0
    Nobs[nj0]=Nobs[nj0]+Nobs0
    


#global limitors, applies for nodes with no observations in adjacent elements 
ZmaxGlobal = +1.
ZminGlobal = -11000.
jMissedMax=np.where(Zmax < ZminGlobal)[0].tolist()
jMissedMin=np.where(Zmin > ZmaxGlobal)[0].tolist()

Zmax[jMissedMax]=ZmaxGlobal
Zmin[jMissedMin]=ZminGlobal

Zmean=np.divide(Zsum , Nobs.astype(float) )

j=np.where(Nobs0==0)[0].tolist()
Zmean[j]= 0. 

flout="MeshNodeInfo."+meshname+flinname+".nc"

with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
    ncout.createDimension('n_nodes' , nn)
    ncout.setncattr("mesh file", mshfl)
    ncout.setncattr("data file", flin)
    ncout.setncattr("imposed global z maximum",ZmaxGlobal)
    ncout.setncattr("imposed global z minimum",ZminGlobal)
    ncout.setncattr("data file", flin)
    ncout.setncattr("SrcFieldType", "unstructured") 
    
    n_var=ncout.createVariable('nobs', 'i8', ('n_nodes',))
    n_var.long_name     = 'number of observations in adjacent elements'
    n_var[:]=Nobs[:]

    zmax_var=ncout.createVariable('zmax', 'f8', ('n_nodes',))
    zmax_var.long_name     = 'maximum observation value in adacent elements'
    zmax_var[:]=Zmax[:]

    zmin_var=ncout.createVariable('zmin', 'f8', ('n_nodes',))
    zmin_var.long_name     = 'minimum observation value in adacent elements'
    zmin_var[:]=Zmin[:]
    
    zmean_var=ncout.createVariable('zmean', 'f8', ('n_nodes',))
    zmean_var.long_name     = 'mean observation value over adacent elements'
    zmean_var[:]=Zmean[:]

    ncout.close
