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

import numpy as np
from scipy.sparse import csr_matrix

def construct_gradient_operators_on_elements(x, y, elements):
    """
    Constructs global df/dx and df/dy operators as CSR sparse matrices 
    evaluating derivatives directly *on* the elements.
    Returns:
    --------
    Dx, Dy : csr_matrix of shape (N_elements, N_nodes)
        Derivative operators mapping node values to element constant derivatives.
    """
    num_nodes = len(x)
    num_elements = len(elements)
    # Each element has 3 nodes -> exactly 3 non-zero entries per row
    # Total entries = num_elements * 3
    I = np.zeros(num_elements * 3, dtype=np.int32)
    J = np.zeros(num_elements * 3, dtype=np.int32)
    val_x = np.zeros(num_elements * 3, dtype=np.float64)
    val_y = np.zeros(num_elements * 3, dtype=np.float64)
    idx = 0
    elements=elements-1
    eind=0
    A=np.zeros(num_elements, dtype=np.float64)
    for e in range(num_elements):
        n_idx = elements[e]
        # Unpack vertex coordinates
        x1, x2, x3 = x[n_idx]
        y1, y2, y3 = y[n_idx]
        # Signed double area via cross product determinant
        two_area =np.abs( (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) )
        A[eind]=two_area/2.
        eind=eind+1
        amin=10.**-10
        if two_area > amin:
        # Standard geometric definitions for P1 linear gradients
        # df/dx = (1 / 2A) * [ (y2-y3)*f1 + (y3-y1)*f2 + (y1-y2)*f3 ]
        # df/dy = (1 / 2A) * [ (x3-x2)*f1 + (x1-x3)*f2 + (x2-x1)*f3 ]
            b1, c1 = y2 - y3, x3 - x2
            b2, c2 = y3 - y1, x1 - x3
            b3, c3 = y1 - y2, x2 - x1
            # Map element row 'e' to the 3 participating global nodes
        else:
            b1,c1,b2,c2,b3,c3 = 0.,0.,0.,0.,0.,0. 
            two_area=1.
        for i, (b, c) in enumerate([(b1, c1), (b2, c2), (b3, c3)]):
            I[idx] = e
            J[idx] = n_idx[i]
            val_x[idx] = b / two_area
            val_y[idx] = c / two_area
            idx += 1
    Dx = csr_matrix((val_x, (I, J)), shape=(num_elements, num_nodes))
    Dy = csr_matrix((val_y, (I, J)), shape=(num_elements, num_nodes))

    return Dx, Dy, A

deg2m=111000.
deg2rad=np.pi/180.
def construct_gradient_operators_on_elements_meters(lon, lat, elements):
    """
    Constructs global df/dx and df/dy operators as CSR sparse matrices 
    evaluating derivatives directly *on* the elements.
    Returns:
    --------
    Dx, Dy : csr_matrix of shape (N_elements, N_nodes)
        Derivative operators mapping node values to element constant derivatives.
    """
    num_nodes = len(lon)
    num_elements = len(elements)
    # Each element has 3 nodes -> exactly 3 non-zero entries per row
    # Total entries = num_elements * 3
    I = np.zeros(num_elements * 3, dtype=np.int32)
    J = np.zeros(num_elements * 3, dtype=np.int32)
    val_x = np.zeros(num_elements * 3, dtype=np.float64)
    val_y = np.zeros(num_elements * 3, dtype=np.float64)
    idx = 0
    elements=elements-1
    eind=0
    A=np.zeros(num_elements, dtype=np.float64)
    for e in range(num_elements):
        n_idx = elements[e]
         # Unpack vertex coordinates
        x1, x2, x3 = lon[n_idx]
        y1, y2, y3 = lat[n_idx]
        latc=(y1+y2+y3)/3.
        # translate to (x1,y1) and convert to meters
        x2=(x2-x1)*np.cos(latc*deg2rad)*deg2m
        x3=(x3-x1)*np.cos(latc*deg2rad)*deg2m
        y2=(y2-y1)*deg2m
        y3=(y3-y1)*deg2m
        x1=0.
        y1=0.
        # Signed double area via cross product determinant
        #two_area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
        two_area = x2*y3 - x3*y2 
        #two_area = np.abs( x1*(y2-y3) + x2*(y3-y1) + x3*(y2-y1) )
        A[eind]=two_area/2.
        eind=eind+1
        amin=10.**-10
        if two_area > amin:
        # Standard geometric definitions for P1 linear gradients
        # df/dx = (1 / 2A) * [ (y2-y3)*f1 + (y3-y1)*f2 + (y1-y2)*f3 ]
        # df/dy = (1 / 2A) * [ (x3-x2)*f1 + (x1-x3)*f2 + (x2-x1)*f3 ]
            b1, c1 = y2 - y3, x3 - x2
            b2, c2 = y3 - y1, x1 - x3
            b3, c3 = y1 - y2, x2 - x1
            # Map element row 'e' to the 3 participating global nodes
        else:
            b1,c1,b2,c2,b3,c3 = 0.,0.,0.,0.,0.,0. 
            two_area=1.
        for i, (b, c) in enumerate([(b1, c1), (b2, c2), (b3, c3)]):
            I[idx] = e
            J[idx] = n_idx[i]
            val_x[idx] = b / two_area
            val_y[idx] = c / two_area
            idx += 1
    Dx = csr_matrix((val_x, (I, J)), shape=(num_elements, num_nodes))
    Dy = csr_matrix((val_y, (I, J)), shape=(num_elements, num_nodes))

    return Dx, Dy, A

#main program
import sys
import numpy as np
import netCDF4 as nc
import xarray as xr
import scipy.sparse as sp

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import cg

nargin = len(sys.argv) - 1
flin=sys.argv[1]
mshfl=sys.argv[2]

x, y, z, e, bnd=FE.loadWW3Mesh(mshfl)

nn=len(x)
flinslash=flin.rfind('/')+1
meshslash=mshfl.rfind('/')+1
meshname=mshfl[meshslash:-3]
flinname=flin[flinslash:-3]

# read in Sparse matrix information
flinSM="SparseMatrix."+meshname+flinname+".nc"
data = nc.Dataset(flinSM,"r")
Nrows = getattr(data, "Nrows")
Ncols = getattr(data, "Ncols")
row=np.array(data["row"][:])
col=np.array(data["col"][:])
S  =np.array(data["S"][:])
# The interpolation matrix from mesh nodes to bathymetry points
A = sp.coo_matrix((S, (row, col-1)), shape=(Nrows, Ncols)).tocsr() # NOTE: col is mesh node number, 1 .. Ncols  

zi = np.array(data["z"][:]) # bathymetry observations
ei = np.array(data["err"][:]) # bathymetry observation standard error

R = diags(1./ei, format="csr") # diagonal matrix of 1 / (obs std err) = sqrt(inv(R)) where R =obs covariance

# construct smoothness penalty terms 
# minimize misfit +  < df/dx ^2 + df/dy ^2 >
# Dx,Dy, EleArea=construct_gradient_operators_on_elements(x, y, e) # gradient operators 
Dx,Dy, EleArea=construct_gradient_operators_on_elements_meters(x, y, e) # gradient operators 
print("EleArea")
print(np.min(EleArea))
e0=e-1
EleArea0=FE.ElementArea(x, y, e0)

W = diags(np.sqrt(EleArea), format="csr") # square root of element area
#W = diags(np.sqrt(EleArea0), format="csr") # square root of element area

WDx=W*Dx
WDy=W*Dy

M=np.transpose(WDx)*WDx + np.transpose(WDy)*WDy # squared slope penalty,spatial integral


# scale penalty with number of observations
#Nobs=len(zi)
#Nnodes=len(x)
#gamma=1.*float(Nobs)/float(Nnodes)

# Use background bathy estimate to set smoothing parameter
z0=-z # background mesh bathymetry
jup=np.where(z0>0.)
z0[jup]=0.
jbad=np.where(np.isnan(z0))
z0[jbad]=0.
r=R*(A*z0-zi)
MP0=np.sum(r**2)
print("MisfitPenalty background = "+str(MP0))
print(M.shape)
print(z0.shape)
SP0=np.dot(z0,M*z0)
print("SmoothPenalty background = "+str(SP0))
#gamma=.05*MP0/SP0 # set smoothing parmater here 
gamma=(10.**-8)*MP0/SP0 # set smoothing parmater here 
print("set smoothing parameter, gamma = "+str(gamma))

AR = R * A
Rz = R * zi
#G= 2.* np.transpose(AR) * AR + 2.*gamma*M
#b =2.* np.transpose(AR)*Rz
G=  np.transpose(AR) * AR + gamma*M
b = np.transpose(AR)*Rz

z0=-z;
jup=np.where(z0>0.)
z0[jup]=0.

x, info = cg(G, b, x0=z0, rtol=10.**(-9))

flroot=meshname+flinname

np.savetxt(flroot+".optbat.txt",x)

output_residual=False
r=A*x-zi
if output_residual:
    np.savetxt(flroot+".residual.txt",r, fmt="%.6e")
SP=np.sqrt(np.mean(r**2))
print("Misfit RMS = "+str(SP) + "(m)")

r=zi-np.mean(zi)
SP=np.sqrt(np.mean(r**2))
print("data RMS  = " + str(SP) + "(m)")

r=R*(A*x-zi)
if output_residual:
    np.savetxt(flroot+".residual_normalized.txt",r, fmt="%.6e")
SP=np.sum(r**2)
print("misfit penalty = "+str(SP))

r=WDx*x
np.savetxt(flroot+".smoothness_penalty.txt",r)
SP=gamma*np.sum(r**2)
print("smoothness penalty = "+str(SP))

################################################################
#Apply local limits to solution:

flinZlim="MeshNodeInfo."+meshname+flinname+".nc"
data = nc.Dataset(flinZlim,"r")
zmean=np.array(data["zmean"][:])
zmin=np.array(data["zmin"][:])
zmax=np.array(data["zmax"][:])

jbad=np.where(zmin>zmax)[0].tolist()
zminGlobal=-11000.
zmaxGlobal= 1.
zmin[jbad]=zminGlobal
zmax[jbad]=zmaxGlobal

G=np.transpose(AR) * AR + gamma*M
b=np.transpose(AR)*Rz

c=-b
from scipy.optimize import minimize, Bounds
# 2. Define the Objective Function, Gradient, and Hessian
# To handle large sparse matrices efficiently, provide analytical grad and hess
def objective(x):
    # 0.5 * x^T * Q * x + c^T * x
    return 0.5 * x.dot(G.dot(x)) + np.dot(c, x)
def gradient(x):
    # Q * x + c
    return G.dot(x) + c
def hessian(x):
    # Returns the sparse CSR matrix directly
    return G

lower_bounds = zmin 
upper_bounds = zmax
bounds = Bounds(lower_bounds, upper_bounds)

# 4. Run Optimization
# Initial guess (e.g., zeros or midpoints)
#x0 = np.zeros(n)
x0=z0
res = minimize(
    fun=objective,
    x0=x0,
    jac=gradient,
    hess=hessian,
    method='trust-constr',
    bounds=bounds,
    options={'maxiter': 10000, 'verbose': 1}
)

# Optimized decision variables
print("Optimization successful:", res.success)
#print("Optimal x:", res.x)
print("Minimum value:", res.fun)

x=res.x

flroot="WithLim."+meshname+flinname

np.savetxt(flroot+".optbat.txt",x)

r=A*x-zi
if output_residual:
    np.savetxt(flroot+".residual.txt",r, fmt="%.6e")
SP=np.sqrt(np.mean(r**2))
print("Misfit RMS = "+str(SP) + "(m)")

r=zi-np.mean(zi)
SP=np.sqrt(np.mean(r**2))
print("data RMS  = " + str(SP) + "(m)")

r=R*(A*x-zi)
if output_residual:
    np.savetxt(flroot+".residual_normalized.txt",r, fmt="%.6e")
SP=np.sum(r**2)
print("misfit penalty = "+str(SP))

r=WDx*x
np.savetxt(flroot+".smoothness_penalty.txt",r)
SP=gamma*np.sum(r**2)
print("smoothness penalty = "+str(SP))

