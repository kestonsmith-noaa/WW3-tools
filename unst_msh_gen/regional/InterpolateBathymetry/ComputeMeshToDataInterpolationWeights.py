import numpy as np
import netCDF4 as nc
import os
import FiniteElementMeshRoutines as FE
from scipy.sparse import diags

deg2rad= np.pi / 180.
deg2kmY=111.
def QuickDistance(lat1, lon1, lats2, lons2):
    deg2kmX=np.cos( deg2rad * lat1 )*deg2kmY
    d = np.min( ((lat1-lats2)*deg2kmY)**2 + ((lon1-lons2)*deg2kmX)**2  )
    d = np.sqrt(d)
    return d

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
                FoundElement[k]=0
                if len(jg)>0:
                   j = jg[0]
                else:
                   j=0
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

def FindElementX(x,y,e,xi,yi,ne):
    j=-9999
    if ne > 1: 
        j=0
        while j < ne:
            xl = np.squeeze(x[e[j,:]])
            yl = np.squeeze(y[e[j,:]])
            IsIn=IsInElement(xl,yl,xi,yi)
            if IsIn :
                return j
            else:
                j=j+1
        if j > ne - 1:
            j=-9999
    else: #one element
        xl = np.squeeze(x[e])
        yl = np.squeeze(y[e])
        IsIn=IsInElement(xl,yl,xi,yi)
        if IsIn:
            j=0
    return j

def compute_mesh_to_point_interp_weights_box(x, y, e, xi, yi):
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

    x = np.asarray(x)
    y = np.asarray(y)
    e = np.asarray(e)
    xi = np.asarray(xi)
    yi = np.asarray(yi)
    
    # convert node indexs to 0 .. nn-1
    e0=e-1
    xc = np.mean(x[e0], axis=1)
    yc = np.mean(y[e0], axis=1)
    
    lleps=10.**-6
    xmin = np.min(x[e0], axis=1)- lleps
    ymin = np.min(y[e0], axis=1)- lleps
    xmax = np.max(x[e0], axis=1)+ lleps
    ymax = np.max(y[e0], axis=1)+ lleps
    # Initialize output arrays
    n_points = len(xi)
    n_elements = len(e)
    weights = np.zeros((n_points, 3)) -99999.
    nodes = np.zeros((n_points, 3), dtype=int) -99999
    elenum = np.zeros(n_points, dtype=int) -99999
    Dist2EleCenter = np.zeros(n_points) -99999.

    FoundElement=np.ones(n_points)
    for k in range(n_points):
        x0=xi[k]
        y0=yi[k]
        jxU=np.where( x0 <= xmax )[0].tolist()
        jxD=np.where( x0 >= xmin )[0].tolist()
        jyU=np.where( y0 <= ymax )[0].tolist()
        jyD=np.where( y0 >= ymin )[0].tolist()
        ji=list( set(jxU) & set(jxD) & set(jyU) & set(jyD) )
        print(str(k)+" : length of searcg elements : " +str(len(ji)) )
        if len(ji) > 0 :
            deg2kmX=np.cos( np.pi * y0 / 180.)*deg2kmY
            distances = np.abs( deg2kmX*(x0 - xc[ji]) + 1j*deg2kmY*(y0 - yc[ji]))
            sindx = np.argsort(distances)
            ji=np.array(ji)
            jg = ji[sindx] # ordered from cloesest to farthest
            j=FindElementX(x,y,e0[jg,:],x0,y0,len(jg))
            if j < 0:
                FoundElement[k]=0
                j = 0 # take cloesest element but note that it is not found
            Dist2EleCenter[k] = distances[j]
            j=jg[j]
            
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
            nodes[k, :] = e[j, :] #<- indexed 1 .. nn

        else :
            FoundElement[k]=0
            j=-99999
        elenum[k] = j # j + 1 ?
        # Progress reporting every 100 iterations
        if (k + 1) % 1000 == 0:
            print(f"Progress: {k+1}/{n_points}")
    return weights, nodes, elenum, Dist2EleCenter, FoundElement


def IsInElementV(x,y,xp,yp):
    npoints=len(xp)
    IsIn=np.zeros(npoints,dtype=int)
    c1 = (x[1] - x[0]) * (yp - y[0]) - (y[1] - y[0]) * (xp - x[0])
    c2 = (x[2] - x[1]) * (yp - y[1]) - (y[2] - y[1]) * (xp - x[1])
    c3 = (x[0] - x[2]) * (yp - y[2]) - (y[0] - y[2]) * (xp - x[2])
    
    j1=np.where( c1 > 0. )[0].tolist()
    j2=np.where( c2 > 0. )[0].tolist()
    j3=np.where( c3 > 0. )[0].tolist()
    jp=list( set(j1) & set(j2) & set(j3) )
    j1=np.where( c1 < 0. )[0].tolist()
    j2=np.where( c2 < 0. )[0].tolist()
    j3=np.where( c3 < 0. )[0].tolist()
    jn=list( set(j1) & set(j2) & set(j3) )
    
    IsIn[jp]=1
    IsIn[jn]=1
    return IsIn

def lenX(obj):
    if hasattr(obj, '__len__'):
        return obj.__len__()
    if isinstance(obj, (int, float, complex)):
        return 1

"""
        deg2kmX=np.cos( np.pi * yc[k] / 180.)*deg2kmY
        distances = np.abs( deg2kmX*(xi - xc[k]) + 1j*deg2kmY*(yi - yc[k]))
        rx= deg2kmX*(np.abs(xl[0]-xl[2] ) + np.abs( xl[2]-xl[1])  + np.abs( xl[1]-xl[0]) )
        ry= deg2kmY*(np.abs(yl[0]-yl[2] ) + np.abs( yl[2]-yl[1])  + np.abs( yl[1]-yl[0]) )
        ji=np.where( distances < rx+ry )
        jf =np.where(FoundElement==0)[0].tolist()
        jxU=np.where( xi[jf] <= xmax[k] )[0].tolist()
        jf=[ jf[jj] for jj in jxU ]
        jxD=np.where( xi[jf] >= xmin[k] )[0].tolist()
        jf=[ jf[jj] for jj in jxD ]
        jyU=np.where( yi[jf] <= ymax[k] )[0].tolist()
        jf=[ jf[jj] for jj in jyU ]
        jyD=np.where( yi[jf] >= ymin[k] )[0].tolist()
        jf=[ jf[jj] for jj in jyD ]
        ji=jf
"""

def compute_mesh_to_point_interp_weights_element_loop(x, y, e, xi, yi):
    
    # Convert to numpy arrays
    from itertools import chain
    deg2kmY=111.

    x = np.asarray(x)
    y = np.asarray(y)
    e = np.asarray(e)
    xi = np.asarray(xi)
    yi = np.asarray(yi)
    
    # convert node indexs to 0 .. nn-1
    e0=e-1
    xc = np.mean(x[e0], axis=1)
    yc = np.mean(y[e0], axis=1)
    
    lleps=10.**-8
    print("lleps="+str(lleps))
    
    xmin = np.min(x[e0], axis=1)- lleps
    ymin = np.min(y[e0], axis=1)- lleps
    xmax = np.max(x[e0], axis=1)+ lleps
    ymax = np.max(y[e0], axis=1)+ lleps
    # Initialize output arrays
    n_points = len(xi)
    n_elements = len(e)
    weights = np.zeros((n_points, 3)) -99999.
    nodes = np.zeros((n_points, 3), dtype=int) -99999
    elenum = np.zeros(n_points, dtype=int) -99999
    Dist2EleCenter = np.zeros(n_points) -99999.

    FoundElement=np.zeros(n_points, dtype=int)
    
    from scipy.spatial import cKDTree
    PC=np.vstack((xc,yc)).T
    PI=np.vstack((xi,yi)).T
    print(PC.shape)
    print(PI.shape)
    tree = cKDTree(PI)
    ne=len(xmin)
    jnf=range(n_points) # j not found
    for k in range(ne):
        xl=x[e0[k,:]]
        yl=y[e0[k,:]]
        rx=np.max( np.abs([xl-xc[k]]) )
        ry=np.max( np.abs([yl-yc[k]]) )
        r_cutoff=np.sqrt(rx**2+ry**2)
        ji = tree.query_ball_point(PC[k,:], r=r_cutoff, workers=-1)
        nji=lenX(ji)
        if nji>0:
            IsIn=IsInElementV(xl,yl,xi[ji],yi[ji])
            npf=sum(IsIn)
            if npf > 0 :
                jfound=np.zeros(npf,dtype=int)
                n=0
                for kk in range(nji):
                    if IsIn[kk]==1:
                        jfound[n]=ji[kk]
                        n=n+1
                FoundElement[jfound]=1
                elenum[jfound] = k
                nodes[jfound,:] = np.ones((npf,1),dtype=np.int32)*e[k,:]
                
                xil=xi[jfound]
                yil=yi[jfound]

                deg2kmX=np.cos( np.pi * yc[k] / 180.)*deg2kmY
                distances = np.abs( deg2kmX*(xil - xc[k]) + 1j*deg2kmY*(yil - yc[k]))
                Dist2EleCenter[jfound]=distances
                
                for j in range(len(jfound)):
                    # Compute barycentric coordinates using shoelace formula for areas
                    # Area of triangle formed by vertices 1, 2, and point (a3)
                    x0=xil[j]
                    y0=yil[j]
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
                    weights[jfound[j], :] = [a1, a2, a3] / total_area
        if k % 1000 == 0:
            print(f"Progress: {k}/{ne}")
            print(str(k)+" : length of search points : " +str(len(ji)) + " of ne="+str(ne) )
            print(str(np.sum(FoundElement))+" : points of " +str(len(FoundElement)) + " found")

    return weights, nodes, elenum, Dist2EleCenter, FoundElement



import os

# write job array to jobcard
def WriteInterpJobscriptSLURM(jobfl, flin, mshfl, Njobs): #, ComputeNodes):
    flinslash=flin.rfind('/')+1
    meshslash=mshfl.rfind('/')+1
    meshname=mshfl[meshslash:len(mshfl)-3]
    TmpOutDir="InterpWeights."+meshname+"."+flin[flinslash:len(flin)-3]
    meshslash=mshfl.rfind('/')+1
    meshname=mshfl[meshslash:len(mshfl)-3]
    with open(jobfl, 'w') as f:
        f.write("#!/bin/bash \n")
        f.write("#SBATCH --job-name=FindInterpWeightsMasterscript \n")
        f.write("#SBATCH --ntasks=1 \n") # ntasks per interpolation
        f.write("#SBATCH --time=08:00:00 \n") 
        f.write("#SBATCH --output=mpi_test_%j.log \n")
        f.write("#SBATCH --error=%j.err \n")
        f.write("#SBATCH --account=marine-cpu \n")
       # f.write("#SBATCH --nodes="+str(ComputeNodes)+" \n")
        f.write("#SBATCH --ntasks-per-core=1"+" \n")
        f.write("#SBATCH --array=0-"+str(Njobs-1)+" \n")
        f.write("#SBATCH --mem=16G \n") # EACH of the 100 tasks gets 4 GB of RAM

        f.write(" \n")
        f.write("source /scratch3/NCEPDEV/climate/Keston.Smith/esmpyenv2/esmpyenv/bin/activate \n")
        f.write("pip list \n")
        current_dir = os.getcwd()
        f.write("cd "+current_dir+"\n")
        f.write("# calculate interpolation weights for data in parallel\n")
        f.write("python -u ComputeMeshToDataInterpolationWeights.py "+flin+" "+mshfl+" $SLURM_ARRAY_TASK_ID " + str(Njobs)+" > InterpJob.$SLURM_ARRAY_TASK_ID.out \n")
        f.write("wait ; \n")
        f.write("echo \"finished\" ; \n")


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
meshslash=mshfl.rfind('/')+1
meshname=mshfl[meshslash:len(mshfl)-3]
TmpOutDir="InterpWeights."+meshname+flin[flinslash:len(flin)-3]

#if nargin==3 then write the parallel jobcard that carries out the interpolation weight calculation
if nargin ==3 :
    Njobs=int(sys.argv[3])
    os.makedirs(TmpOutDir, exist_ok=True)
    try:
        os.remove(TmpOutDir+"/*.txt")
    except:
        print("directory "+TmpOutDir+" is alread empty")
    WriteInterpJobscriptSLURM("jobcardComputeMeshToObsInterpWeightsSLURM",flin,mshfl,Njobs)
    print("Made parallel jobcards to create interpolation weights with "+str(Njobs)+" processes. Next step:")
    print("sbatch jobcardcompute_unstr_to_rwps_interp_weightsSLURM")
    sys.exit()
    
# only continue if nargin == 4
#if nargin==4 then do the part of the mesh for jobID of Njobs
jobID=int(sys.argv[3])
Njobs=int(sys.argv[4])
weights_file = TmpOutDir+"/Part.IntrpWghts."+str(jobID)+".txt"

print("saving output to:"+weights_file)

x, y, z, e, bnd=FE.loadWW3Mesh(mshfl)

nn=len(x)

RWPS_LON=True
if RWPS_LON : 
    j=np.where(x>90.) # should be empty
    x[j]=x[j]-360.

#read point bathymetry data set
data = nc.Dataset(flin,"r")
xi=np.asarray(data["lon"][:])
yi=np.asarray(data["lat"][:])
zi=np.asarray(data["z"][:])
ei=np.array(data["StdErr"][:])

if RWPS_LON :
    j=np.where(xi>90.) # should be empty
    xi[j]=xi[j]-360.

# rudimentary data reduction 
jgm=np.where(  (
            np.logical_and( 
            np.logical_and( xi>=np.min(x), xi<=np.max(x)), 
            np.logical_and( yi>=np.min(y), yi<=np.max(y)) )
            ) )
            
#jg=jg[0].tolist
jgm=jgm[0]
n1=len(jgm)
n0=len(xi)
print("removing data points outside of mesh range. Reduced numnber of data points from "+str(n0)+" to "+str(n1) )
xi=xi[jgm]
yi=yi[jgm]
zi=zi[jgm]
ei=ei[jgm]

#############################################################
# balance node load 
# on PBS system jobID runs from 1 to Njobs
# on SLURM system jobID runs from 0 to Njobs-1 so this partition scheam needs to be modified
nni=len(xi)
xis=np.sort(xi)
js = np.argindices = np.argsort(xi)
mm=round(nni/Njobs)
xil=xis[range(0,nni,mm)] # geographic partitions for parallelization

lonW=xil[jobID]
if jobID==Njobs-1:
    lonE=xis[-1]+1.
else:
    lonE=xil[jobID+1]
if jobID==0:
    lonW=xis[0]-1.

latS=np.min(yi)-1.
latN=np.max(yi)+1.
########################################

e0=e-1
xc = np.mean(x[e0], axis=1)
yc = np.mean(y[e0], axis=1)

nn=len(x)

print(e.shape)
ne=e.shape[0]
print(ne)

xwin=[lonW,lonE]
ywin=[latS,latN]
print("calculating interpolation weights for [W"+str(lonW)+": E "+str(lonE)+": S "+str(latS)+": N "+str(latN)+"]")
#Find target mesh nodes in window
jxUi=np.where( xi < np.max(xwin))[0].tolist()
jxDi=np.where( xi >= np.min(xwin))[0].tolist()
ji=list( set(jxUi) & set(jxDi) )
xi=xi[ji]
yi=yi[ji]
zi=zi[ji]
ei=ei[ji]

point_set=jgm[ji]

#Find source elements in and near window
SearchWidth=1.0 # search rectangle extended beyond window to this number (degrees lat lon)
xmin=np.min(x[e0],axis=1)
xmax=np.max(x[e0],axis=1)
jxU=np.where( xmin < np.max(xwin)+SearchWidth )[0].tolist()
jxD=np.where( xmax > np.min(xwin)-SearchWidth )[0].tolist()
je=list( set(jxU) & set(jxD) )

weights, nodes, elenum, Dist2EleCenter, FoundElement=compute_mesh_to_point_interp_weights_element_loop(x, y, e[je,:], xi, yi)

#recover global element number after localization
jg=np.where(elenum > -1)[0].tolist()
je=np.array(je)
elenum[jg]=je[elenum[jg]]

AllOutMatrix=np.column_stack( ( point_set, zi, ei, nodes, weights, elenum, FoundElement, Dist2EleCenter) )
np.savetxt(weights_file, AllOutMatrix, fmt='%d %f %f %d %d %d %.6f %.6f %.6f %d %d %f')

