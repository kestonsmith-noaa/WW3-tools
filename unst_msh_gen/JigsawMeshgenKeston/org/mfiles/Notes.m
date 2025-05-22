
g=loadmsh('mesh_1kmUSCL_15kmOOFltOcnLLH.msh');
x=g.point.coord(:,1);y=g.point.coord(:,2);z=g.point.coord(:,3);e=g.tria3.index(:,1:3);


ax =

  -82.5444  -50.4751   32.9625   50.4790

  ph=patch_global_local(x,y,z,e',ax);colorbar;