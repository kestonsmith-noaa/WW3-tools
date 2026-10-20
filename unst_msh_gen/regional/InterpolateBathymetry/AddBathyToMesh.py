
#
# Function for taking bathymetry field in a text file and replacing bathymetry
# from a WW3 .msh files with the values from the file specified on the command
# line:
#
# $python3 AddBathyToMesh.py MyBathy.txt 
#
# creates a new WW3 mesh with bathymetry from MyBathy.txt
#
#

import FiniteElementMeshRoutines as FE
import numpy as np
import sys

print("Defining meshname and other variables from file DefineVariables.py")
with open("DefineVariables.py", "r") as file:
    code = file.read()
    print(code)
    exec(code)


meshslash=mesh.rfind('/')+1
mshnm=mesh[meshslash:len(mesh)-4]

flin=sys.argv[1]
flout=flin+".WW3.msh"
nargin = len(sys.argv) - 1
if nargin<2:
    MinBatDepth=float('inf')
else:
    MinBatDepth=float(sys.argv[2])

x, y, z0, e, bnd = FE.loadWW3Mesh(mesh)

nn=len(x)
z=np.zeros(nn)
f=open(flin, 'r')
for k in range(nn):
    line=f.readline()
    z[k]=float(line.strip())
    z[k]=-z[k] # convert to positive is Hdown bathymetric depth
    z[k]=max(z[k],MinBatDepth) 
f.close
FE.WriteWW3Mesh(flout,x,y,z,e,bnd)

