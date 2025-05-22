%from ex_9.m

clear 
close all
flbat='/mnt/sdb/keston/AliMeshTools/RT_GEBCO/Global10km.mat'
flcl='/mnt/sda/keston/GlobalCoastline/gshhg/GSHHS_shp/f/GSHHS_f_L1.shp'
[A,R] = shaperead(flcl);
eval(['load ',flbat]);
bat.z=z;
bat.lon=lon;
bat.lat=lat;

%make pslg;
x=[];y=[];edges=[];
for k=1:length(A)
    minX=min(A(k).X);
    maxX=max(A(k).X);
    minY=min(A(k).Y);
    maxY=max(A(k).Y);
    if and(maxX<0,minX>-100),
        if and(maxY>0,maxY<80),
            Z=A(k).X(1:end-1)+i*A(k).Y(1:end-1);
            if (Z(1)==Z(end))
            x=[x;real(Z(1:end-1))]
        end
    end

end

    initjig;                            % load jigsaw

%------------------------------------ setup files for JIGSAW

    rootpath = fileparts( ...
        mfilename( 'fullpath' )) ;
    rootpath = ...
        fullfile(rootpath, '..') ;

    opts.geom_file = ...                % domain file
        fullfile(rootpath,...
            './', 'proj.msh') ;

    opts.jcfg_file = ...                % config file
        fullfile(rootpath,...
            'cache', 'aust.jig') ;

    opts.mesh_file = ...                % output file
        fullfile(rootpath,...
            'cache', 'mesh.msh') ;

    opts.hfun_file = ...                % sizing file
        fullfile(rootpath,...
            'cache', 'spac.msh') ;

    geom=loadmsh('/home/keston/jigsaw-geo-matlab/files/aust.msh')
 x=geom.point.coord(:,1); y=geom.point.coord(:,2);z=geom.point.coord(:,3);
 edges= geom.edge2.index(:,1:2);
 %chains=edges2chains(edges);
 clf
 plot(x,y,'k.',x(edges'),y(edges'),'r');%slow



    %geom = loadmsh(fullfile( ...
    %    rootpath,'files','aust.msh')) ;

    %topo = loadmsh(fullfile( ...
    %    rootpath,'files','topo.msh')) ;
