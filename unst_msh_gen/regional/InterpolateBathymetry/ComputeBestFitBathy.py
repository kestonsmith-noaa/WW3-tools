import numpy as np
import netCDF4 as nc
import esmpy
import os
HOMErwps = os.environ['HOMErwps']
CIWrwps=HOMErwps+'/dev/compute_interpolation_weights/ush'

def loadWW3Mesh(fl):
    print("mesh file="+fl)
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
        B=A.lstrip()
        values = B.split(" ")
        if len(values)>5:
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
    return xi, yi, ei, zi

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


# Move to the directory where the job was submitted
#-->cd $PBS_O_WORKDIR
#-->echo "Running on node: $(hostname)"
#--> python my_script.py --task $PBS_ARRAY_INDEX

import os
def WriteInterpJobscriptPBS(fl,flin,mshfl,Njobs, ComputeNodes):
    
    meshslash=mshfl.rfind('/')+1
    TmpOutDir="STOFSInterpWeights."+mshfl[meshslash:len(mshfl)-4]
    WghtFl="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".txt"
    WghtFl="InterpWeights."+mshfl[meshslash:len(mshfl)-4]+".stofs.txt"
    WghtFlNetCDF="InterpWeights."+mshfl[meshslash:len(mshfl)-4]+".stofs.nc"
    with open(fl, 'w') as f:
        f.write("#PBS -S /bin/bash\n")
        f.write("#PBS -N UnstrCompIntpWghts.part\n")
        f.write("#PBS -j oe\n")
        f.write("#PBS -q dev\n")
        f.write("#PBS -A NWPS-DEV\n")
        f.write("#PBS -l walltime=08:00:00\n")
#        f.write("#PBS -l select=1:ncpus="+str(Njobs)+":mem=256gb\n")
        f.write("#PBS -l select=1:ncpus=1:mem=32gb\n")
#        f.write("#PBS -l place=excl\n") #makes each sub job exclusive in PBS
        f.write("#PBS -l debug=true\n")
        f.write("#PBS -r y\n")
        f.write("#PBS -J 1-"+str(Njobs)+"\n")

        f.write("module reset\n")
        f.write("module load PrgEnv-intel/8.5.0\n")
        f.write("module load intel/19.1.3.304\n")
        f.write("module load craype/2.7.17\n")
        f.write("module load cray-mpich/8.1.19\n")
        f.write("module load hdf5-C/1.14.0\n")
        f.write("module load netcdf-C/4.9.2\n")
        f.write("module load esmf-C/8.6.0\n")
        f.write("module load ve/hafs/2.1\n")

        f.write("pip list -v\n")

        current_dir = os.getcwd()
        f.write("cd "+current_dir+"\n")

        f.write("# calculate interpolation weights in parallel geographically \n")
        f.write("python "+CIWrwps+"/compute_unstr_to_rwps_interp_weights.py "+flin+" "+mshfl+" $PBS_ARRAY_INDEX " + str(Njobs)+" > InterpJob.$PBS_ARRAY_INDEX.out \n")
#        f.write("python compute_unstr_to_rwps_interp_weights.py "+flin+" "+mshfl+" $PBS_ARRAY_INDEX " + str(Njobs)+" > InterpJob.$PBS_ARRAY_INDEX.out \n")


def WriteInterpJobscriptSLURM(fl,flin,mshfl,Njobs, ComputeNodes):
    
    meshslash=mshfl.rfind('/')+1
    TmpOutDir="STOFSInterpWeights."+mshfl[meshslash:len(mshfl)-4]
    WghtFl="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".txt"
    WghtFlNetCDF="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".nc"
    with open(fl, 'w') as f:
        f.write("#!/bin/bash \n")
        f.write("#SBATCH --job-name=STOFS_interp_masterscript \n")
        f.write("#SBATCH --ntasks=1 \n") # ntasks per interpolation
        f.write("#SBATCH --time=08:00:00 \n") 
        f.write("#SBATCH --output=mpi_test_%j.log \n")
        f.write("#SBATCH --error=%j.err \n")
        f.write("#SBATCH --account=marine-cpu \n")
        f.write("#SBATCH --nodes="+str(ComputeNodes)+" \n")
        f.write("#SBATCH --ntasks-per-core=1"+" \n")
        f.write("#SBATCH --array=0-"+str(Njobs-1)+" \n")

        f.write(" \n")

        f.write("module purge \n")
        f.write("module use /scratch4/NCEPDEV/marine/Ali.Salimi/Hera_Data/HR4-OPT/FromJessica/Keston/ICunstructuredRuns15km-implicit-450s/global-workflow/sorc/ufs_model.fd/modulefiles \n")
        f.write("module load ufs_ursa.intel \n")
        f.write("module load py-scipy/1.14.1 \n")
        f.write("module load py-netcdf4/1.7.1.post2 \n")
        f.write("pip list \n")
        current_dir = os.getcwd()
        f.write("cd "+current_dir+"\n")

        f.write("# calculate interpolation weights in parallel geographically \n")
        f.write("python compute_unstr_to_rwps_interp_weights.py "+flin+" "+mshfl+" $SLURM_ARRAY_TASK_ID " + str(Njobs)+" > InterpJob.$SLURM_ARRAY_TASK_ID.out \n")


#main program
import sys

flin=sys.argv[1]
mshfl=sys.argv[2]

data = nc.Dataset(flin,"r")
xi=np.array(data["lon"][:])
yi=np.array(data["lat"][:])
zi=np.array(data["z"][:])
ei=np.array(data["StdErr"][:])

j=np.where(xi>90.)
xi[j]=xi[j]-360.
j=np.where(zi<0.)
xi=xi[j]
yi=yi[j]
zi=zi[j]

x,y,e,z=loadWW3Mesh(mshfl)
#jx=find(np.where( xi>np.min(x) and xi<np.max(x))
#jy=find(np.where( yi>np.min(y) and yi<np.max(y))
#j=np.where(  ( xi>np.min(x) and xi<np.max(x) )  and ( yi>np.min(y) and yi<np.max(y) ) ) 
#j=np.where(  ( xi>np.min(x) & xi<np.max(x) )  & ( yi>np.min(y) & yi<np.max(y) ) ) 
#j=np.where(  ( (xi>np.min(x), xi<np.max(x) ,  yi>np.min(y),  yi<np.max(y) )  ).all)
bx= np.logical_and( xi>np.min(x), xi<np.max(x) )
by= np.logical_and( yi>np.min(y), yi<np.max(y) )
b= np.logical_and( bx,by )
j=np.where(b)
xi=xi[j]
yi=yi[j]
zi=zi[j]


e0=e-1
e
xc = np.mean(x[e0], axis=1)
yc = np.mean(y[e0], axis=1)
np.savetxt('xc.txt', xc)
np.savetxt('yc.txt', yc)
np.savetxt('xi.txt', xi)
np.savetxt('yi.txt', yi)
np.savetxt('zi.txt', zi)
np.savetxt('ei.txt', ei)

print(np.min(xi))
print(np.max(xi))
print(np.min(x))
print(np.max(x))

weights, nodes, elenum, Dist2EleCenter, FoundElement=compute_mesh_to_mesh_interp_weights(x, y, e, xi, yi)

np.savetxt('weights.txt', weights)
np.savetxt('nodes.txt', nodes, fmt="%d")
np.savetxt('elenum.txt', elenum, fmt="%d" )
np.savetxt('Dist2EleCenter.txt', Dist2EleCenter)
np.savetxt('FoundElement.txt',FoundElement, fmt="%d")
print(weights.shape)



import numpy as np
import netCDF4 as nc
import xarray as xr
import scipy.sparse as sp


weights=np.loadtxt('weights.txt')
nodes=np.loadtxt('nodes.txt',dtype=int)
elenum=np.loadtxt('elenum.txt',dtype=int)
#np.savetxt('Dist2EleCenter.txt', Dist2EleCenter)
FoundElement=np.loadtxt('FoundElement.txt',dtype=int)
print(weights.shape)
#FoundElement=np.int(FoundElement)
print(FoundElement)

jg=np.where(FoundElement==1)
jg
jg=jg[0]
print(jg)

n=0

nd=len(jg)
print(str(nd) + "= nd")
row=np.zeros(nd*3, dtype=int)
col=np.zeros(nd*3, dtype=int)
val=np.zeros(nd*3)

for k in range(nd):
   # print(str(k)+ ":"+str(nd))
    for j in range(3):
        #row[n]=jg[k]
        row[n]=k
        col[n]=nodes[jg[k],j]
        val[n]=weights[jg[k],j]
        n=n+1
n_s=n
        
x,y,e,z=loadWW3Mesh(mshfl)
nn=len(x)

meshslash=mshfl.rfind('/')+1
mshnm=mshfl[meshslash:len(mshfl)-4]
flout=flin[0 : -3]+".weights."+mshnm[0 : -3 ]+".nc"
#flout='IntWghts.nc'

with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
    ncout.createDimension('n_s' , n_s)
    ncout.setncattr("Nrows", nd)
    ncout.setncattr("Ncols", nn)
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

print(row.shape)
print(col.shape)
print(val.shape)

print("max(row)")
print(max(row))
print("max(col)")
print(max(col))

print("min(row)")
print(min(row))
print("min(col)")
print(min(col))

print(nn)
print(nd)
print(n_s)

#matrix = sp.coo_matrix((val, (row-1, col-1)), shape=(nd,nn)).tocsr()
matrix = sp.coo_matrix((val, (row, col-1)), shape=(nd,nn)).tocsr()

zi=np.loadtxt('zi.txt')
zi=zi[jg]
ei=np.loadtxt('ei.txt')
ei=ei[jg]

#zero out land
#jp=np.where(zi>0.)
#zi[jp]=0.
#ignore land
jp=np.where(zi>0.)
ei[jp]=ei[jp]*10000.

print(matrix)
print(zi)
print(zi.shape)
print(x.shape)
print(matrix.shape)

A=matrix

# f(x)=(A x - z )^T (Ax  -z)
# f(x) = x^T A^T A x  -2 (A x)^T z - z^T z
# f(x) = x^T A^T A x  -2 (A x)^T z + c
# f(x) = x^T A^T A x  -2 x A^T z + c
# f(x)     ~= (1/2) x^T G x - b^t x
# where:
# b = 2 * A^Tz
# G = 2 * A^T A

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import cg

G= 2.* np.transpose(A) * A
b =2.* np.transpose(A)*zi
print(G.shape)
print(b.shape)

#x, info = cg(C, b, tol=1e-6)
#x, info = cg(G, b)
x, info = cg(G, b, x0=z,rtol=1e-8)
print(info)
print(x)

np.savetxt(mshnm[0 : -3 ]+".optbat.txt",x)

r=A*x-zi
np.savetxt(mshnm[0 : -3 ]+".residual.txt",r)

diag=range(nd)
R= sp.coo_matrix((1./ei, (diag, diag)), shape=(nd,nd)).tocsr()

# solve
# f(x)=(A x - z )^T COV^(-1) (Ax  -z)
# f(x)=(A x - z )^T diag(1./r/r) (Ax  -z)
# where: r =observation standard err =>ei
# f(x) = [( A x -z ) diag(1./r) ] [( A x -z ) diag(1./r) ]
# f(x) = [( A x -z ) R]^T [( A x -z ) R]   
# f(x) = [ A R x - R z ]^T [ A R x - R z ]
#
# f(x)     ~= (1/2) x^T G x - b^t x
# where:
# b = 2 * A^T R z
# G = 2 * [A R ]^T [A R]

AR= R * A
Rz=R*zi

G= 2.* np.transpose(AR) * AR
b =2.* np.transpose(AR)*Rz

print(G.shape)
print(b.shape)

#x, info = cg(G, b)
x, info = cg(G, b, x0=z,rtol=1e-8)

print(info)
print(x)

np.savetxt(mshnm[0 : -3 ]+".optbatObsErr.txt",x)

r=A*x-zi
np.savetxt(mshnm[0 : -3 ]+".residualObsErr.txt",r)

r=R*(A*x-zi)
np.savetxt(mshnm[0 : -3 ]+".residualObsErrNormalized.txt",r)
