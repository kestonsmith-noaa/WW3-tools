
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
mshnm=mesh[meshslash:len(mesh)-3]

flin=sys.argv[1]
flout=flin+".WW3.msh"

x, y, z0, e, bnd = FE.loadWW3Mesh(mesh)

zmin=1.
nn=len(x)
z=np.zeros(nn)
f=open(flin, 'r')
for k in range(nn):
    line=f.readline()
    #print(str(k)+" "+line)
    z[k]=float(line.strip())
    z[k]=-z[k]
    z[k]=max(z[k],zmin)
f.close
FE.WriteWW3Mesh(flout,x,y,z,e,bnd)

