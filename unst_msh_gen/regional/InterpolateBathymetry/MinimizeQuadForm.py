import numpy as np
import netCDF4 as nc
import os
import FiniteElementMeshRoutines as FE
from scipy.sparse import diags



WidenLocalConstraint=False
Nrelax=1

ComputeLocalLimitSolution = True
ApplyLandValues=True
ImposeLimitsToObs=True
z_coast=-2.0
zminGlobal=-11000.
zmaxGlobal= 1.
zmaxGlobal= z_coast
#Smooth="TPS"
#gamma0=10.**(-4)
Smooth="TVD"
gamma0=10.**(-2)


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

def  construct_stiffness_operator_meters( lon, lat, e ) :

    nn = len(lon)
    ne = len(e)
    deg2m=111000.
    
    idx = 0
    e=e-1 # move index to 0 ... nn-1

    A=np.zeros(ne, dtype=np.float64)
    I=np.zeros(ne*9,dtype=np.int32)
    J=np.zeros(ne*9,dtype=np.int32)
    K=np.zeros(ne*9, dtype=np.float64) # stiffness matrix
    M=np.zeros(ne*9, dtype=np.float64) # mass matrix
    Mlump=np.zeros(nn, dtype=np.float64) # lumped mass (as vector)
    
    dx=np.zeros(3)
    dy=np.zeros(3)
    
    for k in range(ne):
        n_idx = e[k,:]
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

        dx[0]=x2-x3
        dx[1]=x3-x1
        dx[2]=x1-x2
        dy[0]=y2-y3
        dy[1]=y3-y1
        dy[2]=y1-y2
        
        n0=n_idx[0]
        n1=n_idx[1]
        n2=n_idx[2]
        
        A[k]=( x[n0]*dy[0]+x[n1]*dy[1]+x[n2]*dy[2] )/2.
        
        for i in range(3):
            ii=e[k,i]
            Mlump[ii]=Mlump[ii]+A[k]/3.0
            for j in range(3):
                ij=e[k,j]
                I[idx]=ii
                J[idx]=ij
                if ii==ij :
                    M[idx]=A[k]/6.
                else :
                    M[idx]=A[k]/12.

                K[idx]= ( dx[i]*dx[j]+dy[i]*dy[j] ) / A[k] /4.
                idx+=1
                
    K   = csr_matrix((K, (I, J)), shape=(nn, nn))
    Minv = diags( 1./Mlump , format="csr")

    return K, M, Mlump

def construct_gradient_operators_on_nodes_meters(lon, lat, e):

    nn = len(lon)
    ne = len(e)
    deg2m=111000.
    
    idx = 0
    e=e-1 # move index to 0 ... nn-1

    A=np.zeros(ne, dtype=np.float64)
    I=np.zeros(ne*9,dtype=np.int32)
    J=np.zeros(ne*9,dtype=np.int32)
    ddx=np.zeros(ne*9, dtype=np.float64)
    ddy=np.zeros(ne*9, dtype=np.float64)
    K=np.zeros(ne*9, dtype=np.float64) # stiffness matrix
    M=np.zeros(ne*9, dtype=np.float64) # mass matrix
    Mlump=np.zeros(nn, dtype=np.float64)
    
    dx=np.zeros(3)
    dy=np.zeros(3)
    
    for k in range(ne):
        n_idx = e[k,:]
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

        dx[0]=x2-x3
        dx[1]=x3-x1
        dx[2]=x1-x2
        dy[0]=y2-y3
        dy[1]=y3-y1
        dy[2]=y1-y2
        
        n0=n_idx[0]
        n1=n_idx[1]
        n2=n_idx[2]
        
        A[k]=( x[n0]*dy[0]+x[n1]*dy[1]+x[n2]*dy[2] )/2.
        
        for i in range(3):
            ii=e[k,i]
            Mlump[ii]=Mlump[ii]+A[k]/3.0
            for j in range(3):
                ij=e[k,j]
                ddx[idx]=ddx[idx]+dy[j]/6.
                ddy[idx]=ddy[idx]+dx[j]/6.
                I[idx]=ii
                J[idx]=ij
                if ii==ij :
                    M[idx]=A[k]/6.
                else :
                    M[idx]=A[k]/12.

                K[idx]= ( dx[i]*dx[j]+dy[i]*dy[j] ) / A[k] /4.
                idx+=1
                
    DDx = csr_matrix((ddx, (I, J)), shape=(nn, nn))
    DDy = csr_matrix((ddy, (I, J)), shape=(nn, nn))
    K   = csr_matrix((K, (I, J)), shape=(nn, nn))
    M   = csr_matrix((M, (I, J)), shape=(nn, nn))
    Minv = diags( 1./Mlump , format="csr")
    Mlump = diags( Mlump , format="csr")

    return DDx, DDy, M, K, Mlump, Minv

def get_boundary_nodes( elements):
  """Determine boundary nodes of a 2D triangular mesh.

    Parameters:
        nodes: (N, 2) array of (x, y) coordinates for N nodes.
        elements: (M, 3) array of node indices for M triangular elements.

    Returns:
        boundary_nodes: 1D numpy array of unique boundary node indices.
    """
  edge_count = {}

  # Each triangle has 3 edges (pairs of node indices)
  for tri in elements:
    # Ensure consistent edge orientation (sorted tuple) so (i, j) == (j, i)
    edges = [
        tuple(sorted((tri[0], tri[1]))),
        tuple(sorted((tri[1], tri[2]))),
        tuple(sorted((tri[2], tri[0]))),
    ]

    for edge in edges:
      edge_count[edge] = edge_count.get(edge, 0) + 1

  # Boundary edges appear exactly once
  boundary_nodes = set()
  for edge, count in edge_count.items():
    if count == 1:
      boundary_nodes.add(edge[0])
      boundary_nodes.add(edge[1])

  return np.array(sorted(list(boundary_nodes)))

import numpy as np
from scipy.sparse import csr_matrix


def set_rows_to_identity_elements(A, J):
  """Zeros out rows in array J, setting A[j, j] = 1 for each j in J."""
  A = A.tocsr()
  A.sort_indices()

  # Ensure J is a unique NumPy array of integers
  J = np.unique(np.atleast_1d(J).astype(int))

  # 1. Identify slices in A.data corresponding to the rows in J
  # Generate ranges of indices to wipe out
  starts = A.indptr[J]
  ends = A.indptr[J + 1]

  # Vectorized creation of all data/index positions inside target rows
  # equivalent to: np.concatenate([np.arange(s, e) for s, e in zip(starts, ends)])
  counts = ends - starts
  if counts.sum() == 0:
    return A  # No elements to modify in these rows

  # Create flat array of all data coordinates belonging to rows in J
  repeat_offsets = np.repeat(starts, counts)
  sequence_offsets = np.arange(counts.sum()) - np.repeat(
      np.cumsum(counts) - counts, counts
  )
  data_indices = repeat_offsets + sequence_offsets

  # 2. Zero out all existing entries in those rows
  A.data[data_indices] = 0.0

  # 3. Handle setting the diagonal element A[j, j] = 1
  # Find where column index matches row index within the sliced blocks
  col_indices = A.indices[data_indices]
  row_indices_mapped = np.repeat(J, counts)

  diag_mask = col_indices == row_indices_mapped
  actual_diag_data_positions = data_indices[diag_mask]

  if len(actual_diag_data_positions) == len(J):
    # Perfect scenario: All A[j,j] slots already existed structurally
    A.data[actual_diag_data_positions] = 1.0
  else:
    # Safe fallback if some diagonal items structurally didn't exist yet:
    # Convert to LIL or reconstruct cleanly
    A_lil = A.tolil()
    A_lil[J, :] = 0
    A_lil[J, J] = 1
    A = A_lil.tocsr()

  A.eliminate_zeros()
  return A


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

if ImposeLimitsToObs:
    jh=np.where(zi>zmaxGlobal)
    jl=np.where(zi<zminGlobal)
    zi[jh]=zmaxGlobal
    zi[jl]=zminGlobal

R = diags(1./ei, format="csr") # diagonal matrix of 1 / (obs std err) = sqrt(inv(R)) where R =obs covariance

# construct smoothness penalty terms 
# TPS smoothness penalty < (d2f/dx2 )^2 +2*(d2f/dxdy )^2 (d2f/dy2 )^2 +  >
# where d2f/dx2 is the second derivative of f in the x direction, etc. 
match Smooth:
    case "TPS":
#    DDx, DDy, M, K, Mlump, Minv = construct_gradient_operators_on_nodes_meters( x, y, e )
        K, Mass, Mlump = construct_stiffness_operator_meters( x, y, e )
        Minv = diags(1./Mlump, format="csr")
        M  = K * Minv * K
        Msqrt=diags(np.sqrt(1./Mlump), format="csr")*K # M= Msqrt^T * Msqrt
        Dx, Dy, EleArea=construct_gradient_operators_on_elements_meters(x, y, e) # gradient operators 
        W = diags(np.sqrt(EleArea), format="csr") # square root of element area
        WDx=W*Dx
        WDy=W*Dy
    case "TVD":
        Dx, Dy, EleArea=construct_gradient_operators_on_elements_meters(x, y, e) # gradient operators 
        print("EleArea")
        print(np.min(EleArea))
        W = diags(np.sqrt(EleArea), format="csr") # square root of element area
        WDx=W*Dx
        WDy=W*Dy
        M=np.transpose(WDx)*WDx + np.transpose(WDy)*WDy # squared slope penalty,spatial integral
        Msqrt=WDx+WDy

print("M")
print(M)

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
#gamma=(10.**-4)*MP0/SP0 # set smoothing parmater here 
gamma=gamma0*MP0/SP0 # set smoothing parmater here 
print("set smoothing parameter, gamma = "+str(gamma))


jl=np.where(zi > z_coast)
zi[jl]=z_coast

AR = R * A
Rz = R * zi

G=  np.transpose(AR) * AR + gamma*M
b = np.transpose(AR)*Rz


if ApplyLandValues :
    # simple determine land boundaries
    bnd=get_boundary_nodes( e )
    print("bnd")
    print(bnd)
    bnd=bnd-1 # indexed 1...nn -> 0 ...(nn-1)
    r=.1
    jE=np.where( x[bnd] > np.min(x)+r )[0].tolist()
    jW=np.where( x[bnd] < np.max(x)-r )[0].tolist()
    jN=np.where( y[bnd] > np.min(y)+r )[0].tolist()
    jS=np.where( y[bnd] < np.max(y)-r )[0].tolist()
    jint=list( set(jE) & set(jW) & set(jN) & set(jS) )
    bnd=bnd[jint]
    G=set_rows_to_identity_elements(G, bnd)
    b[bnd]=z_coast

z0=-z;
jup=np.where(z0>0.)
z0[jup]=0.

x, info = cg(G, b, x0=z0, rtol=10.**(-8))
print(info)
print("info")

import osqp
from scipy import sparse
from scipy.sparse import csr_array, csc_array

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

match Smooth:
    case "TVD":
        rx=WDx*x
        ry=WDy*x
        SP=gamma*np.sum(rx**2+ry**2)
        r=np.sqrt(rx**2+ry**2) # I think this is right
    case "TPS":
        r=Msqrt*x

np.savetxt(flroot+".smoothness_penalty.txt",r)
SP=gamma*np.sum(r**2)
print("smoothness penalty = "+str(SP))

################################################################
#Apply local limits to solution:
ComputeLocalLimitSolution = True
if ComputeLocalLimitSolution :
    flinZlim="MeshNodeInfo."+meshname+flinname+".nc"
    data = nc.Dataset(flinZlim,"r")
    zmean=np.array(data["zmean"][:])
    zmin=np.array(data["zmin"][:])
    zmax=np.array(data["zmax"][:])

    jbad=np.where(zmin>zmax)[0].tolist()
    zmin[jbad]=zminGlobal
    zmax[jbad]=zmaxGlobal



    if WidenLocalConstraint:
        zl=np.vstack((zmin,zmax))
        np.savetxt('zlimits0.txt',zl.T)
        ne=len(e)
        print("ne = " + str(ne))
        for jit in range(Nrelax):
            zminE=np.zeros(ne)+zmaxGlobal
            zmaxE=np.zeros(ne)+zminGlobal
            for k in range(ne):
                indx=e[k,:]-1
                zmaxE[k]=np.max(zmax[indx])
                zminE[k]=np.min(zmin[indx])
            zmin0=zmin
            zmax0=zmax
            for k in range(ne):
                indx=e[k,:]-1
                for ii in range(3):
                    zmin0[indx[ii]]=min(zmin0[indx[ii]],zminE[k])
                    zmax0[indx[ii]]=max(zmax0[indx[ii]],zmaxE[k])
            zmin=zmin0
            zmax=zmax0
        
        zl=np.vstack((zmin,zmax))
        zlE=np.vstack((zminE,zmaxE))
        np.savetxt('zlimits.txt',zl.T)
        np.savetxt('zlimitsE.txt',zlE.T)






    G=np.transpose(AR) * AR + gamma*M
    b=np.transpose(AR)*Rz

    print("a")
    from scipy.optimize import lsq_linear
    if ApplyLandValues :
        # simple determine land boundaries
        bnd=get_boundary_nodes( e )
        bnd=bnd-1 # indexed 1...nn -> 0 ...(nn-1)
        print("bnd")
        print(bnd)
        print(len(bnd))
        r=.1 # assume oriented rectangular boundary
        jE=np.where( x[bnd] > np.min(x)+r )[0].tolist()
        jW=np.where( x[bnd] < np.max(x)-r )[0].tolist()
        jN=np.where( y[bnd] > np.min(y)+r )[0].tolist()
        jS=np.where( y[bnd] < np.max(y)-r )[0].tolist()
        jint=list( set(jE) & set(jW) & set(jN) & set(jS) )
        bnd=bnd[jint]
        print(bnd)
        print(len(bnd))
        G=set_rows_to_identity_elements(G, bnd)
        b[bnd]=z_coast
        zmax[bnd]=z_coast+1.
        zmin[bnd]=z_coast-1.
        # this may not be correct- need to remove land points from the
        # matrix and account on RHS 


    # 4. Run Optimization
    # Initial guess (e.g., zeros or midpoints)
    #x0 = np.zeros(n)
    x0=z0
    import osqp
    from scipy import sparse
    from scipy.sparse import csr_array, csc_array
    prob = osqp.OSQP()
    Aeye = diags(np.ones(nn), format="csr")
# wrong!!!    prob.setup(G.tocsc(), b, Aeye.tocsc(), zmin, zmax, alpha=1.0, eps_abs=1e-7, eps_rel=1e-7 )
#    prob.setup(G.tocsc(), -b, Aeye.tocsc(), zmin, zmax, alpha=1.0, eps_abs=1e-7, eps_rel=1e-7 )
    prob.setup(G.tocsc(), -b, Aeye.tocsc(), zmin, zmax, alpha=1.0, eps_abs=1e-7, eps_rel=1e-7 , max_iter=10000)
    
    #x_cg, info = cg(G, b, x0=z0, rtol=10.**(-8))
    #prob.warm_start(x=x_cg) #converge more slowly to same solution
    #prob.warm_start(x=z0) 
    
    res = prob.solve()
    x=res.x

    flroot="WithLim."+meshname+flinname

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

    match Smooth:
        case "TVD":
            rx=WDx*x
            ry=WDy*x
            SP=gamma*np.sum(rx**2+ry**2)
            r=np.sqrt(rx**2+ry**2) # I think this is right
        case "TPS":
            r=Msqrt*x

    np.savetxt(flroot+".smoothness_penalty.txt",r)
    SP=gamma*np.sum(r**2)
    print("smoothness penalty = "+str(SP))


