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
meshslash=mshfl.rfind('/')+1
meshname=mshfl[meshslash:len(mshfl)-3]

TmpOutDir="InterpWeights."+meshname+flin[flinslash:len(flin)-3]
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
#    if len(jg) > 0:
    AllData=np.vstack((AllData,AllOutMatrix[jg,:]))

AllData=AllData[1 : ,:]
#np.savetxt("AllData.txt", AllData, fmt='%d %f %f %d %d %d %.6f %.6f %.6f %d %d %f')
#365042 -5528.000000 10.000000 3812 6516 8606 0.360925 0.471654 0.167421 7607 1 17.212554

flout="InterpWeights."+meshname+flin[flinslash:len(flin)-3]+".nc"

N,twelve=AllData.shape
#                                   0       1    2    [3,4,5] [6,7,8]   9       10               11
#AllOutMatrix=np.column_stack( ( point_set, zi,  ei,   nodes, weights, elenum, FoundElement, Dist2EleCenter) )
#np.savetxt(weights_file, AllOutMatrix, fmt='%d %f %f %d %d %d %.6f %.6f %.6f %d %d %f')
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

flout="SparseMatrix."+meshname+flin[flinslash:len(flin)-3]+".nc"

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

"""    
    c_var=ncout.createVariable('col', 'i4', ('n_s',))
    c_var.long_name     = 'column index'
    c_var[:]=col[:]
    
    s_var=ncout.createVariable('S', 'f8', ('n_s',))
    s_var.long_name     = 'matrix value'
    s_var[:]=val[:]
"""

"""    
    print(data_point_number)
    new_points=AllOutMatrix[jg,0]
    print(new_points)
    new_points=new_points.astype(int)
    print(data_point_number)    
    print(new_points)
    new_points=new_points.tolist
    print(new_points)
    data_point_number.extend(new_points)
"""    
"""
#                                   0       1    2    [3,4,5] [6,7,8]   9       10               11
#AllOutMatrix=np.column_stack( ( point_set, zi,  ei,   nodes, weights, elenum, FoundElement, Dist2EleCenter) )
#np.savetxt(weights_file, AllOutMatrix, fmt='%d %f %f %d %d %d %.6f %.6f %.6f %d %d %f')
408182 -4582.000000 10.000000 5805 5808 5804 0.163625 0.613954 0.222421 11620 1 70.732705
408183 -3441.000000 10.000000 5805 5808 5804 0.259573 0.615523 0.124904 11620 1 76.439576
408184 -4537.000000 10.000000 5805 5808 5804 0.355521 0.617092 0.027388 11620 1 92.139210
408185 -4586.000000 10.000000 5806 5808 5805 0.075990 0.550091 0.373919 11624 1 71.793484

flout=TextWeightFl[0:len(TextWeightFl)-4]+".nc"

xi, yi, ei, zi = iutil.loadWW3Mesh(mshfl)
nn_dst=len(xi)
data = nc.Dataset(flin,"r")
#read spaital dimensions and determine if input mesh is curvilinear or regular
x=np.asarray(data["x"][:])
#convert to RWPS coordinates

jEast=np.where(x>90)
x[jEast]=x[jEast]-360.

nn_src=len(x)

#Read in weights from cat of text output(rows are not in order)
F=np.loadtxt(TextWeightFl)

node_dst=F[:,0] # row number (destination node)
ele_num_src=F[:,1] # Source element number (not used)
node_src=F[:,2:5] # Source node
DistToCenterKM=F[:,5] # distance of destination node from source element center (not used) 
weights=F[:,6:9] # interpolation weights from source nodes to destination node

node_dst=node_dst.astype(int)
node_src=node_src.astype(int)
ele_num_src=ele_num_src.astype(int)

nnz=F.shape[0]
n_s=3*nnz
row=np.zeros(n_s,dtype=int)
col=np.zeros(n_s,dtype=int)
val=np.zeros(n_s)

n=0
for k in range(nnz):
    for j in range(3):
        row[n]=node_dst[k]
        col[n]=node_src[k,j]
        val[n]=weights[k,j]
        n=n+1
print("n="+str(n))
print("n_s="+str(n_s))
n_s=n

with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
    ncout.createDimension('n_s' , n_s)
    ncout.setncattr("Nrows", nn_dst)
    ncout.setncattr("Ncols", nn_src)
    ncout.setncattr("SrcFieldType", "unstructured") 
    
    r_var=ncout.createVariable('row', 'i4', ('n_s',))
    r_var.long_name     = 'row index'
    r_var[:]=row[:]
    
    c_var=ncout.createVariable('col', 'i4', ('n_s',))
    c_var.long_name     = 'column index'
    c_var[:]=col[:]
    
    s_var=ncout.createVariable('S', 'f8', ('n_s',))
    s_var.long_name     = 'matrix value'
    s_var[:]=val[:]

#Consider adding (x,y) destination and (x,y) for source. This is usefull for extrapolation
    if AddExtrapolationSupport:
        x=np.asarray(data["x"][:])
        y=np.asarray(data["y"][:])
        ncout.createDimension('nn_src' , nn_src)
        ncout.createDimension('nn_dst' , nn_dst)
        
        xsrc_var=ncout.createVariable('x_src', 'f8', ('nn_src',))
        xsrc_var.long_name     = 'interpolation source node longitude'
        xsrc_var[:]=x[:]
        
        ysrc_var=ncout.createVariable('y_src', 'f8', ('nn_src',))
        ysrc_var.long_name     = 'interpolation source node latitude'
        ysrc_var[:]=y[:]
        
        xdst_var=ncout.createVariable('x_dst', 'f8', ('nn_dst',))
        xdst_var.long_name     = 'interpolation destination node longitude'
        xdst_var[:]=xi[:]
    
        ydst_var=ncout.createVariable('y_dst', 'f8', ('nn_dst',))
        ydst_var.long_name     = 'interpolation destination node latitude'
        ydst_var[:]=yi[:]

        zdst_var=ncout.createVariable('depth_dst', 'f8', ('nn_dst',))
        zdst_var.long_name     = 'interpolation destination node latitude'
        zdst_var[:]=zi[:]

Extrapolate=sys.argv[4]

if Extrapolate:
    dist2bnd=np.full(len(xi), np.inf) #all points are inside boundary- No boundary with this type of extrapolation
else:
    matrix = sp.coo_matrix((weights, (row-1, col-1)), shape=(nn_src,nn_dst)).tocsr()
    row_sum = matrix.sum(axis=1)
    j0=np.where( row_sum==0 ) # destination nodes with no coverage from interpolation matrix
    j0=np.array(j0[0]).tolist()
    u0=np.ones(xi.shape)
    nan=float("nan")
    u0[j0]=nan
    dist2bnd=iutil.CalculateDistanceToInterpEnvelope(xi,yi,u0, 1.)

dist2bnd_file = "DistToBndy."+mshfl[meshslash:len(mshfl)-4]+".stofs.nc"
with nc.Dataset(dist2bnd_file, 'w', format='NETCDF4') as ncout:
    ncout.createDimension('node' , nn_dst)
    d_var=ncout.createVariable('dist2bnd', 'f4', ('node',))
    d_var.long_name     = 'distance to boundary'
    d_var.units         = 'km'
    d_var.standard_name = 'distance to boundary'
    d_var[:]=dist2bnd[:]
    
    z_var=ncout.createVariable('depth', 'f4', ('node',))
    z_var.long_name     = 'mesh depth'
    z_var.units         = 'm'
    z_var.standard_name = 'depth'
    z_var[:]=zi[:]

"""
