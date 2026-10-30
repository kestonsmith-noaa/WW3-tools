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
ComputeWeights=True
if nargin>2 :
    if int(sys.argv[3])==0:
        ComputeWeights=False

weights=np.loadtxt('weights.txt')
nodes=np.loadtxt('nodes.txt',dtype=int)
elenum=np.loadtxt('elenum.txt',dtype=int)
#np.savetxt('Dist2EleCenter.txt', Dist2EleCenter)
FoundElement=np.loadtxt('FoundElement.txt',dtype=int)
zi=np.loadtxt('zi.txt')
ei=np.loadtxt('ei.txt')

print("weights.shape")
print(weights.shape)
print(FoundElement)

print("data zmax = " +str(np.max(zi)) )
print("data zmin = " +str(np.min(zi)) )
print("data zmean = " +str(np.mean(zi)) )
ExcludeLand=False
if ExcludeLand :
    jneg=np.where( zi <= 0)
    jg  =np.where( FoundElement==1 )
    jg=jg[0]
    #jg=list(set(jg) & set(jneg))
    print("jg")
    print(jg)
    jg=np.array(np.intersect1d(jg,jneg))
    print(jg)
#    jg=np.where( np.logical_and( zi < 0., FoundElement==1 ) )
else:
    jg=np.where(FoundElement==1)
    jg=jg[0]


nd=len(jg)
print(str(nd) + "= nd")

print("jg")
print(jg)

n=0
row=np.zeros(nd*3, dtype=int)
col=np.zeros(nd*3, dtype=int)
val=np.zeros(nd*3)
for k in range(nd):
    for j in range(3):
        row[n]=k
        col[n]=nodes[jg[k],j]
        val[n]=weights[jg[k],j]
        n=n+1
n_s=n

#zi=np.loadtxt('zi.txt')
zi=zi[jg]
#ei=np.loadtxt('ei.txt')
ei=ei[jg]
print("data zmax filtered= " +str(np.max(zi)) )

x, y, z, e, bnd=FE.loadWW3Mesh(mshfl)

#x,y,e,z=loadWW3Mesh(mshfl)
nn=len(x)

meshslash=mshfl.rfind('/')+1
mshnm=mshfl[meshslash:len(mshfl)-4]
flout=flin[0 : -3]+".weights."+mshnm[0 : -3 ]+".nc"

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


matrix = sp.coo_matrix((val, (row, col-1)), shape=(nd,nn)).tocsr()


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

# minimize
# f(x)=(A x - z )^T (Ax  -z)
# f(x) = x^T A^T A x  -2 (A x)^T z - z^T z
# f(x) = x^T A^T A x  -2 (A x)^T z + c
# f(x) = x^T A^T A x  -2 x A^T z + c
# f(x)     ~= (1/2) x^T G x - b^t x
# where:
# b = 2 * A^Tz
# G = 2 * A^T A

G= 2.* np.transpose(A) * A
b =2.* np.transpose(A)*zi
print(G.shape)
print(b.shape)
#Set initial guess at bathy based on current mesh z
z0=-z;
jup=np.where(z0>0.)
z0[jup]=0.
tol=10^-6
#x, info = cg(C, b, tol=1e-6)
x, info = cg(G, b)
#x, info = cg(G, b, x0=z0,rtol=tol)
print(info)
print(x)

np.savetxt(mshnm+".optbat.txt",x)

r=A*x-zi
np.savetxt(mshnm+".residual.txt",r)

diag=range(nd)
R= sp.coo_matrix((1./ei, (diag, diag)), shape=(nd,nd)).tocsr()

# minimize
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

#x, info = cg(G, b)
z0=-z;
jup=np.where(z0>0.)
z0[jup]=0.
x, info = cg(G, b, x0=z0)
#x, info = cg(G, b, x0=z0,rtol=tol)
np.savetxt(mshnm+".optbatObsErr.txt",x)

r=A*x-zi
np.savetxt(mshnm+".residualObsErr.txt",r)

r=R*(A*x-zi)
np.savetxt(mshnm+".residualObsErrNormalized.txt",r)

# introduce reguralization: penalize df/dx ^2 + df/dy ^2
x, y, z, e, bnd=FE.loadWW3Mesh(mshfl)

Dx,Dy, EleArea=construct_gradient_operators_on_elements(x, y, e)

np.savetxt('EleArea.txt',EleArea)

cooy = Dy.tocoo()
Dy3 = np.column_stack((cooy.row, cooy.col, cooy.data))
np.savetxt('Dy.txt',Dy3)

coox = Dx.tocoo()
Dx3 = np.column_stack((coox.row, coox.col, coox.data))
np.savetxt('Dx.txt',Dx3)

M=np.transpose(Dx)*Dx + np.transpose(Dy)*Dy # squared slope penalty

# Note the other variant is to minimize < df/dx ^2 + df/dy ^2 >
# in which case we want:
# M=np.transpose(Dx)*diag(EleArea)*Dx + np.transpose(Dy)*diag(EleArea)Dy # squared slope penalty

gamma=10.**-6

# minimize rms_misfit + d/dx ^2 + d/dy ^2 
# g(x)=(A f - z )^T COV^(-1) (A f  -z) + lamda (f^t Dx^T Dx  f + f^t Dy^T Dy  f)
# so: (f^t Dx^t Dx  f + f^t Dy^t Dy  f)= f^t (Dx^t Dx + Dy^t Dy) f = f^t M f
# let: M=(Dx^t Dx + Dy^t Dy)
# g(x) = [ A R f - R z ]^t [ A R f - R z ] + lambda f^t M f
#      = c+ f^t [(R^t A^t A R] f - 2 * [A R x]^t z + lambda f^t M f
#      = c+ f^t [(R [A^t A] R + lambda*M] f - 2 [A^T R z]^t f   
# g(x)     ~= (1/2) x^T G x - b^t x
# where:
# b = 2 * A^t R z
# G = 2 * [A R ]^t [A R] + 2 * lambda*M


G= 2.* np.transpose(AR) * AR + 2.*gamma*M
b =2.* np.transpose(AR)*Rz

#x, info = cg(G, b)
z0=-z;
jup=np.where(z0>0.)
z0[jup]=0.
x, info = cg(G, b, x0=z0)
#x, info = cg(G, b, x0=z0,rtol=tol)
np.savetxt(mshnm+".optbatReg.txt",x)
r=A*x-zi
np.savetxt(mshnm+".residualReg.txt",r)
r=R*(A*x-zi)
np.savetxt(mshnm+".residualRegNormalized.txt",r)

# minimize < df/dx ^2 + df/dy ^2 >
ne, three=e.shape
diag=range(ne)
e0=e-1
EleArea0=FE.ElementArea(x, y, e0)

W = diags(np.sqrt(EleArea), format="csr")
WDx=W*Dx
WDy=W*Dy

#M=np.transpose(Dx)*W*Dx + np.transpose(Dy)*W*Dy # squared slope penalty,spatial integral
M=np.transpose(WDx)*WDx + np.transpose(WDy)*WDy # squared slope penalty,spatial integral
print(M)

gamma=10.**(-9)
gamma=10.**(-2) # shows effect
gamma=10.**(-3) # shows effect
gamma=10.**(-2) # shows effect

G= 2.* np.transpose(AR) * AR + 2.*gamma*M
b =2.* np.transpose(AR)*Rz

z0=-z;
jup=np.where(z0>0.)
z0[jup]=0.
#x, info = cg(G, b, x0=z0,rtol=tol)
x, info = cg(G, b, x0=z0)

np.savetxt(mshnm+".optbatRegI.txt",x)

r=A*x-zi
np.savetxt(mshnm+".residualRegI.txt",r)

SP=np.sum(r**2)
print("Misfit RMS^2 = "+str(SP))

r=R*(A*x-zi)
np.savetxt(mshnm+".residualRegINormalized.txt",r)

SP=np.sum(r**2)
print("MisfitPenalty = "+str(SP))

r=WDx*x
np.savetxt(mshnm+".SmoothPenalty.txt",r)
SP=gamma*np.sum(r**2)
print("SmoothPenalty = "+str(SP))
