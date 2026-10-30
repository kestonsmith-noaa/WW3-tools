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
        #two_area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
        two_area =np.abs( (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) )
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
meshname=mshfl[meshslash:len(mshfl)-3]
flout=flin[0 : -3]+".weights."+meshname[0 : -3 ]+".nc"

# read in Sparse matrix information
flin="SparseMatrix."+meshname+flin[flinslash:len(flin)-3]+".nc"
data = nc.Dataset(flin,"r")
Nrows = getattr(data, "Nrows")
Ncols = getattr(data, "Ncols")
row=np.array(data["row"][:])
col=np.array(data["col"][:])
S  =np.array(data["S"][:])
# The interpolation matrix from mesh nodes to bathymetry points
A = sp.coo_matrix((S, (row, col-1)), shape=(Nrows, Ncols)).tocsr() # NOTE: col is mesh node number, 1 .. Ncols  

zi = np.array(data["z"][:]) # bathymetry observations
ei = np.array(data["err"][:]) # bathymetry observation standard error

#diag=range(len(ei))
#R = sp.coo_matrix((1./ei, (diag, diag)), shape=(nd,nd)).tocsr()
R = diags(ei, format="csr") # diagonal matrux of obs std err

# construct smoothness penalty terms 
# minimize misfit +  < df/dx ^2 + df/dy ^2 >
Dx,Dy, EleArea=construct_gradient_operators_on_elements(x, y, e) # gradient operators 

e0=e-1
EleArea0=FE.ElementArea(x, y, e0)

W = diags(np.sqrt(EleArea), format="csr") # square root of element area
W = diags(np.sqrt(EleArea0), format="csr") # square root of element area

WDx=W*Dx
WDy=W*Dy

M=np.transpose(WDx)*WDx + np.transpose(WDy)*WDy # squared slope penalty,spatial integral

gamma=10.**(-9)
gamma=10.**(-2) # shows effect
gamma=10.**(-3) # shows effect
gamma=10.**(-4) # shows effect
#gamma=10.**(-2) # shows effect

AR = R * A
Rz = R * zi
G= 2.* np.transpose(AR) * AR + 2.*gamma*M
b =2.* np.transpose(AR)*Rz

z0=-z;
jup=np.where(z0>0.)
z0[jup]=0.
#x, info = cg(G, b, x0=z0,rtol=tol)
x, info = cg(G, b, x0=z0)
np.savetxt(meshname+"optbat.txt",x)

r=A*x-zi
np.savetxt(meshname+"residual.txt",r)
SP=np.sum(r**2)
print("Misfit RMS^2 = "+str(SP))

r=R*(A*x-zi)
np.savetxt(meshname+"residualNormalized.txt",r)
SP=np.sum(r**2)
print("MisfitPenalty = "+str(SP))

r=WDx*x
np.savetxt(meshname+"SmoothPenalty.txt",r)
SP=gamma*np.sum(r**2)
print("SmoothPenalty = "+str(SP))
