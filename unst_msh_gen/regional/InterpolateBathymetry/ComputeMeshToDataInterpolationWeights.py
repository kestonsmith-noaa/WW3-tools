import numpy as np
import netCDF4 as nc
import os
import FiniteElementMeshRoutines as FE
from scipy.sparse import diags

import concurrent.futures

def loadWW3Mesh(fl):
    f=open(fl, 'r')
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of nodes
    nn=int(header)
    print("nn = "+str(nn))
    xi=np.zeros(nn)
    yi=np.zeros(nn)
    zi=np.zeros(nn)
    k=0
    for i in range(nn):
        A = f.readline()
        values = A.split(" ")
        if len(values)>4:
            xi[k]=values[2]
            yi[k]=values[4]
            zi[k]=values[6]
        else:
            xi[k]=values[1]
            yi[k]=values[2]
            zi[k]=values[3]
        k=k+1
    print("number of nodes read: "+str(k))
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of elements
    ne=int(header)#includes boundary nodes and actual elements
    print("ne="+str(ne)+" -includes boundary nodes")
    nbnd=0
    bnd=[]
    eix=np.zeros((ne,3), dtype=int)
    k=0
    for i in range(ne):
        A = f.readline()
        values = A.split(" ")
        if len(values) == 6:
            if int(values[2])==2:
                bnd.append(int(values[5]))
                nbnd=nbnd+1
        if len(values)>15:
            eix[k,0]=int(values[12])
            eix[k,1]=int(values[14])
            eix[k,2]=int(values[16])
            k=k+1
        elif len(values)>7:
            eix[k,0]=int(values[6])
            eix[k,1]=int(values[7])
            eix[k,2]=int(values[8])
            k=k+1
    ei=eix[range(k),:]
    print("number of open boundary nodes read: "+str(nbnd))
    print("number of elements read: "+str(k))
    return xi, yi, zi, ei, bnd


deg2rad= np.pi / 180.
deg2kmY=111.
def QuickDistance(lat1, lon1, lats2, lons2):
    deg2kmX=np.cos( deg2rad * lat1 )*deg2kmY
    d = np.min( ((lat1-lats2)*deg2kmY)**2 + ((lon1-lons2)*deg2kmX)**2  )
    d = np.sqrt(d)
    return d

def WriteInterpolationWeightsToNetCDF(weights_file,row,col,weights,Nrows,Ncols):
    #create a esmpy style sparse matrix netcdf file
    n_s=len(weights)

    if isinstance(row , list):
        row      = np.array( row  )
    if isinstance( col, list):
        col      = np.array( col )
    if isinstance( weights, list):
        weights = np.array( weights )

    with nc.Dataset(weights_file, 'w', format='NETCDF4') as ncout:
            
        ncout.createDimension('n_s' , n_s)
        ncout.setncattr("Nrows", Nrows)
        ncout.setncattr("Ncols", Ncols)
            
        r_var=ncout.createVariable('row', 'i4', ('n_s',))
        r_var.long_name     = 'row index'
        r_var[:]=row[:]
            
        c_var=ncout.createVariable('col', 'i4', ('n_s',))
        c_var.long_name     = 'column index'
        c_var[:]=col[:]
            
        s_var=ncout.createVariable('S', 'f4', ('n_s',))
        s_var.long_name     = 'matrix value'
        s_var[:]=weights[:]

def IsInElement(x,y,xp,yp):
    IsIn=False
    c1 = (x[1] - x[0]) * (yp - y[0]) - (y[1] - y[0]) * (xp - x[0])
    c2 = (x[2] - x[1]) * (yp - y[1]) - (y[2] - y[1]) * (xp - x[1])
    c3 = (x[0] - x[2]) * (yp - y[2]) - (y[0] - y[2]) * (xp - x[2])
    if ( ( c1 > 0 and c2 > 0 and c3 > 0) or ( c1 < 0 and c2 < 0 and c3 < 0) ):
        IsIn=True
    if (c1*c2*c3 == 0): # include points on triangle
        IsIn=True
    return IsIn

def FindElement(x,y,e,xi,yi):
    nd=len(e.shape)
    e=np.squeeze(e)
    j=-9999
    ne=e.shape[0]
    if nd > 1: 
        xc = np.squeeze(np.mean(x[e], axis=1))
        yc = np.squeeze(np.mean(y[e], axis=1))
        DistanceToElements = np.abs((xi + 1j*yi) - (xc + 1j*yc))
        n=0
        while (j < 0 and n < ne)  :
            n=n+1
            j = np.argmin( DistanceToElements )
            xl = np.squeeze(x[e[j,:]])
            yl = np.squeeze(y[e[j,:]])
            IsIn=IsInElement(xl,yl,xi,yi)
            if IsIn :
                return j
            else:
                DistanceToElements[j]=float('inf')
                j=-9999
    else:
        xl = np.squeeze(x[e])
        yl = np.squeeze(y[e])
        IsIn=IsInElement(xl,yl,xi,yi)
        if IsIn:
            j=0
            return j
        else:
            j=-9999
            return j
        
    return j


def compute_mesh_to_mesh_interp_weights(x, y, e, xi, yi):
    """
    Compute mesh-to-mesh interpolation weights using barycentric coordinates.
    
    Parameters:
    -----------
    x : array-like
        X coordinates of source mesh nodes
    y : array-like
        Y coordinates of source mesh nodes
    e : array-like
        Element connectivity matrix (n_elements x 3 for triangular elements)
        (nodes indexed from 1 to nn)
    xi : array-like
        X coordinates of target points
    yi : array-like
        Y coordinates of target points
        
    Returns:
    --------
    weights : ndarray
        Interpolation weights (n_points x 3)
    nodes : ndarray
        Node indices for each point (n_points x 3)
    elenum : ndarray
        Element number for each point (n_points,)
    Dist2EleCenter : ndarray
        Distance to element center for each element
    """
    
    # Convert to numpy arrays
    
    deg2kmY=111.
    UseNearestEle = False #if True just use closest element center
    N = 12 # search N nearest elements (nearest by element center to target node distance)

    x = np.asarray(x)
    y = np.asarray(y)
    e = np.asarray(e)
    xi = np.asarray(xi)
    yi = np.asarray(yi)

    # convert node indexs to 0 .. nn-1
    e0=e-1
    e
    print(e0)
    xc = np.mean(x[e0], axis=1)
    yc = np.mean(y[e0], axis=1)
    
    # Initialize output arrays
    n_points = len(xi)
    n_elements = len(e)
    weights = np.zeros((n_points, 3))
    nodes = np.zeros((n_points, 3), dtype=int)
    elenum = np.zeros(n_points, dtype=int)
    Dist2EleCenter = np.zeros(n_points)
    
    #UseNearestEle=True
    FoundElement=np.ones(n_points)
    for k in range(n_points):
        x0=xi[k]
        y0=yi[k]
        deg2kmX=np.cos( np.pi * y0 / 180.)*deg2kmY
        distances = np.abs( deg2kmX*(x0 - xc) + 1j*deg2kmY*(y0 - yc))
        if UseNearestEle:
            jg=np.argmin(distances)
            j=jg
        else:
            raw_indices = np.argpartition(distances, N)[:N]
            sorted_sub_indices = np.argsort(distances[raw_indices])
            jg = raw_indices[sorted_sub_indices]
            j=FindElement(x,y,e0[jg,:],x0,y0)
            if j < 0:
                print("couldn't find element for point : "+ str(k))
                print(str(x0)+" : "+str(y0))
                FoundElement[k]=0
                if len(jg)>0:
                   j = jg[0]
                else:
                   j=0
                print("using closest element at distance : "+ str(distances[j]))
            else:
                j=jg[j]

        elenum[k] = j
        Dist2EleCenter[k] = distances[j]
        xl = x[e0[j, :]]
        yl = y[e0[j, :]]
            
        # Compute barycentric coordinates using shoelace formula for areas
        # Area of triangle formed by vertices 1, 2, and point (a3)
        xt = np.array([xl[0], xl[1], x0, xl[0]])
        yt = np.array([yl[0], yl[1], y0, yl[0]])
        a3 = -np.dot(xt[1:4] - xt[0:3], yt[0:3] + yt[1:4]) / 2
        # Area of triangle formed by vertices 3, 1, and point (a2)
        xt = np.array([xl[2], xl[0], x0, xl[2]])
        yt = np.array([yl[2], yl[0], y0, yl[2]])
        a2 = -np.dot(xt[1:4] - xt[0:3], yt[0:3] + yt[1:4]) / 2
        # Area of triangle formed by point, vertices 2, 3 (a1)
        xt = np.array([x0, xl[1], xl[2], x0])
        yt = np.array([y0, yl[1], yl[2], y0])
        a1 = -np.dot(xt[1:4] - xt[0:3], yt[0:3] + yt[1:4]) / 2

        # Normalize to get barycentric weights
        total_area = a1 + a2 + a3
        weights[k, :] = [a1, a2, a3] / total_area

        # Handle extrapolation cases where target is outside element
        # with nearest element center extrapolation
        if np.min(weights[k, :])<0:
            weights[k, :]=1/3

        nodes[k, :] = e[j, :] #<- indexed 1 .. nn
        # Progress reporting every 100 iterations
        if (k + 1) % 1000 == 0:
            print(f"Progress: {k+1}/{n_points}")
    return weights, nodes, elenum, Dist2EleCenter, FoundElement

import os

def WriteInterpJobscriptSLURM(jobfl,flin,mshfl,Njobs)#, ComputeNodes):
    
    WghtFl="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".txt"
    WghtFlNetCDF="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".nc"
    
    flinslash=flin.rfind('/')+1
    TmpOutDir="InterpWeights."+flin[flinslash:len(flin)-3]
    
    weights_file = TmpOutDir+"/Part.IntrpWghts."+str(jobID)+".txt"
    
    with open(jobfl, 'w') as f:
        f.write("#!/bin/bash \n")
        f.write("#SBATCH --job-name=FindInterpWeightsMasterscript \n")
        f.write("#SBATCH --ntasks=1 \n") # ntasks per interpolation
        f.write("#SBATCH --time=08:00:00 \n") 
        f.write("#SBATCH --output=mpi_test_%j.log \n")
        f.write("#SBATCH --error=%j.err \n")
        f.write("#SBATCH --account=marine-cpu \n")
        #f.write("#SBATCH --nodes="+str(ComputeNodes)+" \n")
        f.write("#SBATCH --ntasks-per-core=1"+" \n")
        f.write("#SBATCH --array=0-"+str(Njobs-1)+" \n")
        f.write(" \n")
        f.write("source /scratch3/NCEPDEV/climate/Keston.Smith/esmpyenv2/esmpyenv/bin/activate \n")
        f.write("pip list \n")
        
        current_dir = os.getcwd()
        
        f.write("cd "+current_dir+"\n")
        f.write("# calculate interpolation weights for data in parallel\n")
        f.write("python ComputeMeshToDataInterpolationWeights.py "+flin+" "+mshfl+" $SLURM_ARRAY_TASK_ID " + str(Njobs)+" > InterpJob.$SLURM_ARRAY_TASK_ID.out \n")
        
        JointFile="InterpWeights."+meshname+"."+flin[flinslash:len(flin)-3]+".txt"
        FullWeightsFile=JointFile[0 : -3]
        
        f.write("wait; #wait for parallel array jobs to finish \n")
        f.write("cat "+TmpOutDir+"/Part.IntrpWghts.*.txt > "+ JointFile + "\n")
        f.write("python convert_weights_to_netcdf.py "+JointFile+" "+FullWeightsFile+" \n+")
        f.write("submit minimization script with "+FullWeightsFile)

#main program
import sys
import numpy as np
import netCDF4 as nc
import xarray as xr
import scipy.sparse as sp

from scipy.sparse import csr_matrix

nargin = len(sys.argv) - 1
flin=sys.argv[1]
mshfl=sys.argv[2]

flinslash=flin.rfind('/')+1
TmpOutDir="InterpWeights."+flin[flinslash:len(flin)-3]

#if nargin==3 then write the parallel jobcard that carries out the interpolation weight calculation
if nargin ==3 :
    Njobs=int(sys.argv[3])
    os.makedirs(TmpOutDir, exist_ok=True)
    try:
        os.remove(TmpOutDir+"/*.txt")
    except:
        print("directory "+TmpOutDir+" is alread empty")
    WriteInterpJobscriptSLURM("jobcardcompute_unstr_to_rwps_interp_weightsSLURM",flin,mshfl,Njobs) #, 1)
    print("Made parallel jobcards to create interpolation weights with "+str(Njobs)+" processes. Next step:")
    print("sbatch jobcardcompute_unstr_to_rwps_interp_weightsSLURM")
    sys.exit()

#if nargin==4 then do the part of the mesh for jobID of Njobs
if nargin == 4:
    jobID=int(sys.argv[3])
    Njobs=int(sys.argv[4])
    weights_file = TmpOutDir+"/Part.IntrpWghts."+str(jobID)+".txt"

print("saving output to:"+weights_file)

x, y, e, z = iutil.loadWW3Mesh(mshfl)

nn=len(x)
if RWPS_LON : 
    j=np.where(x>90.) # should be empty
    x[j]=x[j]-360.

#read point bathymetry data set
data = nc.Dataset(flin,"r")
xi=np.asarray(data0["lon"][:])
yi=np.asarray(data0["lat"][:])
zi=np.asarray(data0["z"][:])
ei=np.array(data["StdErr"][:])

if RWPS_LON :
    j=np.where(xi>90.) # should be empty
    xi[j]=x[j]-360.

nni=len(xi)

# divide point bathymetry data set for jobID=0 ... Njobs-1
def split_with_numpy(n, m):
  return [list(arr) for arr in np.array_split(range(n), m)]
all_job_point_sets=split_with_numpy(nni, Njobs)
point_set=all_job_point_sets[jobID]
xi=xi[point_set]
yi=yi[point_set]
zi=zi[point_set]
ei=ei[point_set]

if nargin == 4:
    weights, nodes, elenum, Dist2EleCenter, FoundElement=compute_mesh_to_mesh_interp_weights(x, y, e, xi, yi)
    AllOutMatrix=np.column_stack( (nodes, weights, point_set, zi, FoundElement, Dist2EleCenter) )
    np.savetxt(weights_file, AllOutMatrix.T, fmt='%d %d %d %.6f %.6f %.6f %d %f %d %f')

