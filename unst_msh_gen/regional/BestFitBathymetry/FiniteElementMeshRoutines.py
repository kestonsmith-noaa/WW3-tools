#
# Functions for reading and writing WW3 meshes, calculating element 
# area, and estimateing mesh length scale.
#
# x and y are length nn node coordinates
# e is a (Number-of-Elements by 3) element array
# bnd is a length nb list of open ocean boundary nodes
#   Note:(nodes are assumed to be indexed from 1,2,...,nn within e
#         and bnd)
#
# loadWW3MeshCoords(fl): read mesh nodes and elements from file
#
# loadWW3Mesh(fl):  read mesh nodes, elements, and open boundary nodes
#                   from WW3 .msh file
#
# WriteWW3Mesh(flo,x,y,z,e,bnd): write WW3 .msh file
#
# lengthscale(x, y, e): compute length scale for each element
#
# ElementArea(x, y, e): compute area for each element
#
# ComputeNodeLengthScale(x,y,e): estimate length scale at each node
#

import numpy as np
import math
#from geopy import distance
import re

# Simple distance approximation in kilometers for calculating element
# area and lengthscale.
deg2km=111.132954
deg2rad=np.pi/180.
def Distance(x0,y0,x1,y1):
    d=np.sqrt( 
        ( deg2km*np.cos(deg2rad*(y0+y1)/2 )*(x1-x0) )**2 + 
        ( deg2km*(y1-y0))**2  ) 
    return d

def loadWW3MeshCoords(fl):
    f=open(fl, 'r')
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of nodes
    nn = re.findall(r'\d+', header)
    nn=int(nn[0])
    print(nn)
    xi=np.zeros(nn)
    yi=np.zeros(nn)
    k=0
    for i in range(nn):
        A = f.readline()
        values = A.split(" ")
        if len(values)>4:
            xi[k]=values[2]
            yi[k]=values[4]
        else:
            xi[k]=values[1]
            yi[k]=values[2]
        k=k+1
    print("number of nodes read: "+str(k))
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of elements
    ne=int(header)
    print("ne=",str(ne))
    ei=np.zeros((ne,3), dtype=int)
    k=0
    for i in range(ne):
        A = f.readline()
        values = A.split(" ")
        #need to factor bnd nodes
        if len(values)>15:
            ei[k,0]=int(values[12])
            ei[k,1]=int(values[14])
            ei[k,2]=int(values[16])
            k=k+1
        elif len(values)>7:
            ei[k,0]=int(values[6])
            ei[k,1]=int(values[7])
            ei[k,2]=int(values[8])
            k=k+1
    #ei=ei[0:k,:]    
    print("number of elements read: "+str(k))
    return xi, yi, ei


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
#    ei=eix[range(k),:]
    ei=eix[range(k-1),:]
    print("number of open boundary nodes read: "+str(nbnd))
    print("number of elements read: "+str(k))
    return xi, yi, zi, ei, bnd

def WriteWW3Mesh(flo,x,y,z,e,bnd):
    f=open(flo, 'w')
    f.write("$MeshFormat\n")
    f.write("2 0 8\n")
    f.write("$EndMeshFormat\n")
    f.write("$Nodes\n")
    nn=len(x)
    f.write(str(nn)+"\n")
    for k in range(nn):
        f.write(f"{k+1:8d} {x[k]:6f} {y[k]:6f} {z[k]:5f} \n")
    ne=e.shape[0]
    nbnd=len(bnd)
    f.write("$EndNodes\n")
    f.write("$Elements\n")
    f.write(str(ne+nbnd)+"\n")
    for k in range(nbnd):
        f.write(str(k+1)+" 15 2 0 0 " + str(bnd[k])+"\n")
    for k in range(ne):
        f.write(str(nbnd+k+1)+" 2 3 0 "+str(k+1)+" 0 "+str(e[k,0])+" "+str(e[k,1])+" "+str(e[k,2])+"\n")
    f.write("$EndElements\n")
    f.close

# compute length scale for each element e
def lengthscale(x, y, e):
    ne=e.shape[0]
    lengthscaleE=np.zeros(ne)
    for k in range(ne):
        if k % 10000 ==0:
            print(k)
        x1=x[e[k,0]-1]
        y1=y[e[k,0]-1]
        x2=x[e[k,1]-1]
        y2=y[e[k,1]-1]
        x3=x[e[k,2]-1]
        y3=y[e[k,2]-1]
        D3=Distance(x1,y1,x2,y2)
        D1=Distance(x2,y2,x3,y3)
        D2=Distance(x3,y3,x1,y1)
        lengthscaleE[k]=(D1+D2+D3)/3  
        # mean edge length
    return lengthscaleE

# compute area for each element e
def ElementArea(x, y, e):
    ne=e.shape[0]
    AreaE=np.zeros(ne)
    for k in range(ne):
        if k % 10000 ==0:
            print(k)
        x1=x[e[k,0]-1];
        y1=y[e[k,0]-1];
        x2=x[e[k,1]-1];
        y2=y[e[k,1]-1];
        x3=x[e[k,2]-1];
        y3=y[e[k,2]-1];
        D3=Distance(x1,y1,x2,y2)
        D1=Distance(x2,y2,x3,y3)
        D2=Distance(x3,y3,x1,y1)
        S=(D1+D2+D3)/2
        A=np.sqrt(S * (S - D1) * (S - D2) * (S - D3))
        AreaE[k]=A
    return AreaE

# compute area for each element e
def ElementAreaP(x, y, e):
    ne=e.shape[0]
    AreaE=np.zeros(ne)
    for k in range(ne):
        if k % 10000 ==0:
            print(k)
        x1=x[e[k,0]-1];
        y1=y[e[k,0]-1];
        x2=x[e[k,1]-1];
        y2=y[e[k,1]-1];
        x3=x[e[k,2]-1];
        y3=y[e[k,2]-1];
        D3=np.sqrt( (x1-x2)**2 + (y1-y2)**2 )
        D1=np.sqrt( (x2-x3)**2 + (y2-y3)**2 )
        D2=np.sqrt( (x3-x1)**2 + (y3-y1)**2 )
        S=(D1+D2+D3)/2
        A=np.sqrt(S * (S - D1) * (S - D2) * (S - D3))
        AreaE[k]=A
    return AreaE

# compute area for each element e
def GradientOnElements(x, y, f, e):
    A=ElementAreaP(x, y, e)
    print("length area="+str(len(A)))
    ne=e.shape[0]
    print("ne="+str(ne))
    dfdx=np.zeros(ne)
    dfdy=np.zeros(ne)
    for k in range(ne):
        if k % 10000 ==0:
            print(k)
        i1=e[k,0]-1
        i2=e[k,1]-1
        i3=e[k,2]-1
        dfdx[k]= ( f[i1]*(y[i2]-y[i3])+f[i2]*(y[i3]-y[i1])+f[i3]*(y[i1]-y[i2]) )/(2.*A[k])
        dfdy[k]=-( f[i1]*(x[i2]-x[i3])+f[i2]*(x[i3]-x[i1])+f[i3]*(x[i1]-x[i2]) )/(2.*A[k])
    return dfdx,dfdy

#A=csr_matrix((Avalues, (row_ind, col_ind-1)), shape=(n, nn)).toarray()

from scipy.sparse import csr_matrix

def GradientMatrixToEle(x, y, e):
# Create CSR matrix to compute gradients on elements
    A=ElementAreaP(x, y, e)
    ne=e.shape[0]
    nn=len(x)
    #n=np.max(e)+1
    ddxV=np.zeros(ne*3)
    ddyV=np.zeros(ne*3)
    i=np.zeros(ne*3)
    j=np.zeros(ne*3)
    for k in range(ne):
        if k % 10000 ==0:
            print(k)
            
        j1=e[k,0]-1
        j2=e[k,1]-1
        j3=e[k,2]-1
        
        i[3*k+0]=k
        i[3*k+1]=k
        i[3*k+2]=k
        
        j[3*k+0]=j1
        j[3*k+1]=j2
        j[3*k+2]=j3
        
        ddxV[3*k+0]=+(y[j2]-y[j3])/(2.*A[k])
        ddxV[3*k+1]=+(y[j3]-y[j1])/(2.*A[k])
        ddxV[3*k+1]=+(y[j1]-y[j2])/(2.*A[k])
        
        ddyV[3*k+0]=-(x[j2]-x[j3])/(2.*A[k])
        ddyV[3*k+1]=-(x[j3]-x[j1])/(2.*A[k])
        ddyV[3*k+1]=-(x[j1]-x[j2])/(2.*A[k])
    print("ne="+str(ne))
    print("nn="+str(nn))
    print("i")
    print(np.min(i))
    print(np.max(i))
    print("j")
    print(np.min(j))
    print(np.max(j))
    print("e0")
    print(np.min(e[:,0]))
    print(np.max(e[:,0]))
    print("e1")
    print(np.min(e[:,1]))
    print(np.max(e[:,1]))
    print("e2")
    print(np.min(e[:,2]))
    print(np.max(e[:,2]))
    print("min j="+str(np.min(j) ))
    ddx=csr_matrix((ddxV, (i, j)), shape=(ne, nn)).toarray()
    ddy=csr_matrix((ddyV, (i, j)), shape=(ne, nn)).toarray()
#raise ValueError(f'axis {i} index {idx.max()} exceeds ValueError: axis 1 index 450 exceeds matrix dimension 450

    return ddx,ddy

        
from scipy.sparse import diags

def GradientMatrixToNode(x, y, e):
# Create CSR matrix to compute gradients on elements
    ddxE,ddyE=GradientMatrixToEle(x, y, e)
    nn=len(x)
    ne=e.shape[0]
    A=ElementAreaP(x, y, e)
    
    ddx=csr_matrix(([],([],[])),shape=(nn, nn)).toarray()
    ddy=csr_matrix(([],([],[])),shape=(nn, nn)).toarray()
    W=np.zeros(nn)

    for k in range(ne):
        if k % 10000 ==0:
            print(k)
        j1=e[k,0]-1
        j2=e[k,1]-1
        j3=e[k,2]-1
        for i in range(3):
            j=e[k,i]-1
            ddx[j,j1]=ddx[j,j1]+ddxE[k,j1]*A[k]/3.
            ddx[j,j2]=ddx[j,j2]+ddxE[k,j2]*A[k]/3.
            ddx[j,j3]=ddx[j,j3]+ddxE[k,j3]*A[k]/3.
            ddy[j,j1]=ddy[j,j1]+ddyE[k,j1]*A[k]/3.
            ddy[j,j2]=ddy[j,j2]+ddyE[k,j2]*A[k]/3.
            ddy[j,j3]=ddy[j,j3]+ddyE[k,j3]*A[k]/3.
            W[j]=W[j]+A[k]/3.
        
    WD = diags(1./W)
    WDcsr=WD.tocsr()
    ddx=WDcsr @ ddx
    ddy=WDcsr @ ddy
    

    return ddx,ddy



# approximate mesh lengthscale at mesh nodes.
# lengthscale is the weighted average (by area) of surrounding elements
def ComputeNodeLengthScale(x,y,e):
    nn=np.max(e[:])
    print("ComputeNodeLengthScale nn="+str(nn))
    areaE=ElementArea(x, y, e)
    LengthScaleE=lengthscale(x, y, e)
    
    ne=e.shape[0]
    areaT=np.zeros(nn)
    lsN=np.zeros(nn)
    for k in range(ne):
        for j in range(3):
            i=e[k,j]-1
            lsN[i]=lsN[i]+LengthScaleE[k]*areaE[k]
            areaT[i]=areaT[i]+areaE[k]
    for k in range(nn):
        lsN[k]=lsN[k]/areaT[k]
    return lsN

