mesh = loadmsh('mesh_1kmUSCL_15kmOO.msh');
load  GlobalCoastline5km.mat
topo=BoxSmoothTopo('RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc',3);
meshP=GlobalMesh2Planer(mesh,topo)
meshO = CutOutLand(meshP,S,3)
savemsh('mesh_1kmUSCL_15kmOO_Ocean.msh',meshO);

