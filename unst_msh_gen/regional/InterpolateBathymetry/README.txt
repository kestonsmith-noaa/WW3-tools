These are subroutines to interpolate bathymetry data to a unstructured
meshes with triangular elements.  The interpolation uses Gauss Markov 
smoothing with variable lengthscale proportional to the mesh lengthscale.
The interpolation is computed in parallel across a user determined 
number of mesh divisions.  The program uses an arbitrary number of 
bathymetry data files.  It was designed to use a low resolution file 
with global coverage as well as local high resolution gridded bathymetry,
and global scattered high resolution bathymetry data, such as global
Satellite Derived Bathymetry (SDB). Bathymetry data sets need to be 
downloaded and preprocessed to use a common format. The software is not 
optimized for speed.

The software was motivated by the development of the Regional Wave Prediction System (RWPS).
Here some important regions require blending SDB with low resolution global database.  In regions
where high resolution gridded bathymetry is available, this is prioritized, however data points 
from all databases are included if they are near the point interpolated to. This is intended
to prevent discontinuities at the boundary of data sets. 


To interpolate bathymetry to a mesh:
    ( -1 ) set python environment to include libraries: numpy  sys os 
           argparse time netCDF4 math itertools subprocess scipy and re
           on ursa: $source /scratch3/NCEPDEV/climate/Keston.Smith/esmpyenv2/esmpyenv/bin/activate
           should work 

   *(0) Download bathymetry data to be used (notes are in  
        MakeCommonNetCDFFileFormat.py) 
       
   *(1) python MakeCommonNetCDFFileFormat.py 

   *(2) python SmoothSubsampleNetCDFFiles.py 

   *(3) python SmoothSubsampleCMEMS.py

   *(Skip Steps 0-3 if running on Ursa - Path to existing, processed, bathymetry data sets is set in DefineVariables.py
     The bathymetry data set collection used for RWPS is fairly large ~100 GB )

    (4) python3 DivideMeshNodes.py NP
            NP=number of interpolation jobs to run in parallel

    (5) sbatch jobcardInterpBathy2Mesh
            Runs interpolation job in parallel in slurm environmant

        This step produces ascii files with an entry for each mesh node :    

HawaiiTest.GMU.txt : Bathymetry (m) at nodes estimated using Gauss-Markov smoothing with unknown mean, i.e. ordinary kriging.
HawaiiTest.GMU.std.txt : Posterior standard error (m) for bathymetry estimate in HawaiiTest.GMU.txt
HawaiiTest.GMN.txt : Bathymetry (m) at nodes estimated using Gauss-Markov with nearest data point used for local mean (simple kriging).
                     This estimator seems to perform better around narrow atolls and other small scale shallow features surrounded by deep 
                     waters- at least with the data sets used here. 
HawaiiTest.GMN.std.txt : Posterior standard error (m) for bathymetry estimate in HawaiiTest.GMN.txt
HawaiiTest.InvDist.txt : Bathymetry (m) estimate using inverse distance weighting (provided for error checking)
HawaiiTest.ClosestValue.txt : Bathymetry (m) for nearest data point to node (provided for error checking)
HawaiiTest.NDataPoints.txt : Number of bathymetry data points used in estimate at each node
HawaiiTest.LengthScale.txt : Interpolation length scale (m) at each node- porportional to mesh length scale

    (6) python AddBathyToMesh.py HawaiiTest.GMU.txt
            creates HawaiiTest.GMU.txt.WW3.msh, a WW3 .msh mesh format
            Bathymetry HawaiiTest.GMU.txt is based on Gauss-Markov smoothing with unknown mean, i.e. ordinary kriging.
            Alternatinly you can use HawaiiTest.GMN.txt: Gauss-Markov with nearest data point used for mean
            or 

HawaiiTest.GMU.txt : Bathymetry (m) at nodes estimated using Gauss-Markov smoothing with unknown mean, i.e. ordinary kriging.
HawaiiTest.GMU.std.txt : Posterior standard error (m) for estimate in HawaiiTest.GMU.txt
HawaiiTest.GMN.txt : Bathymetry (m) at nodes estimated using Gauss-Markov with nearest data point used for local mean (simple kriging )
                     Note: This estimate avoids some artifacts that may occur in HawaiiTest.GMU.txt in regions of sharp transitions.
HawaiiTest.GMN.std.txt : Posterior standard error (m) for estimate in HawaiiTest.GMN.txt 
HawaiiTest.InvDist.txt : Bathymetry (m) estimate using inverse distance weighting (provided for error checking)
HawaiiTest.NDataPoints.txt : Number of bathymetry data points used in estimate at each node
HawaiiTest.ClosestValue.txt : Bathymetry (m) for nearest data point to node 
HawaiiTest.LengthScale.txt : Interpolation length scale (m) at each node, porportional to mesh length scale

to change the mesh interpolated to, bathymetry files used, or paramaters values,
edit file DefineVariables.py .

DefineVariables.py
AddBathyToMesh.py
FiniteElementMeshRoutines.py
GaussMarkov.py
InterpolateCRM.partscat.py
KnitOutputBackTogether.py
MakeCommonNetCDFFileFormat.py
SmoothSubsampleCMEMS.py
SmoothSubsampleNetCDFFiles.py

