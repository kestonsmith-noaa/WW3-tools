

import numpy as np
import netCDF4 as nc
import xarray as xr
import scipy.sparse as sp

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
#np=len(xi)
#nn=len(x)
#np=len(elenum)

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
        
mshfl="meshes/rwps.oc_20km_300km.msh"
x,y,e,z=loadWW3Mesh(mshfl)
nn=len(x)


#flout=flin[0 : -3]+".weights."+mshfl[0 : -3 ]+".nc"
flout='IntWghts.nc'

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
ei[jp]=ei[jp]*100000.

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
x, info = cg(G, b)
print(info)
print(x)

np.savetxt("optbat.txt",x)

r=A*x-zi
np.savetxt("residual.txt",r)


diag=range(nd)
R= sp.coo_matrix((1./ei, (diag, diag)), shape=(nd,nd)).tocsr()

#AR= np.transpose(A) * R
AR= R * A
Rz=R*zi


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


G= 2.* np.transpose(AR) * AR
b =2.* np.transpose(AR)*Rz

print(G.shape)
print(b.shape)

#x, info = cg(C, b, tol=1e-6)
x, info = cg(G, b)

print(info)
print(x)

np.savetxt("optbatObsErr.txt",x)

r=A*x-zi
np.savetxt("residualObsErr.txt",r)

