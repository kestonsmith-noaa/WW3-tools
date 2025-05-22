function GenateGlobalWUScoast
close all ; initjig ;

% DEMO-2 -- generate a regionally-refined global grid with a
%   high-resolution "patch" (@37.5KM) embedded within a
%   uniform background grid (@150.KM).
%   JIGSAW is combined with a bisection procedure to improve
%   the regularity of the grid topology.

%------------------------------------ setup files for JIGSAW
rootpath='./'
  %  rootpath = fileparts( ...
  %      mfilename( 'fullpath' )) ;

    opts.geom_file = ...                % domain file
        fullfile(rootpath,...
            'cache','geom.msh') ;

    opts.jcfg_file = ...                % config file
        fullfile('./',...
            'cache','opts.jig') ;

    opts.mesh_file = ...                % output file
        fullfile(rootpath,...
            'cache','mesh.msh') ;

    opts.hfun_file = ...                % sizing file
        fullfile(rootpath,...
            'cache','spac.msh') ;

%------------------------------------ define JIGSAW geometry

    geom.mshID = 'ELLIPSOID-MESH' ;
    geom.radii = 6371 * ones(3,1) ;

   % opts.geom_file='geom'
    savemsh (opts.geom_file,geom) ;


%------------------------------------ compute HFUN over GEOM

%    topo = loadmsh(fullfile( ...
%        rootpath, 'files', 'topo.msh' ) );
    topoJS = loadmsh('./topo.msh');
    clear topo
    fl='RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc'
    lon=ncread(fl,'lon');
    lat=ncread(fl,'lat');
    z=ncread(fl,'bed_elevation');

       w=ones(5,5)/25;
       [nx,ny] = size(z)
       z5=conv2(z,w,'same');
       z5=z5(3:5:nx,3:5:ny);
       lon5=lon(3:5:nx);
       lat5=lat(3:5:ny);

   % topo.point.coord{:,1}=lon(1:end-1);
   % topo.point.coord{:,2}=lat(1:end-1);%fix to centers later
   % topo.value=double(z');
  
    topo.point.coord{:,1}=lon5;
    topo.point.coord{:,2}=lat5;%fix to centers later

    topo.mshID='ELLIPSOID-GRID'
    topo.fileV=3;
    topo.value=double(z5');
    
    xpos = topo.point.coord{1};
    ypos = topo.point.coord{2};
    zlev = reshape( topo.value,length(ypos), length(xpos));

   [XPOS,YPOS] = meshgrid(xpos,ypos);

    XPOS = XPOS * pi/180 ;
    YPOS = YPOS * pi/180 ;

    S = shaperead('us_coastline/tl_2023_us_coastline.shp');
    xt=[];
    yt=[];
    for k=1:length(S);
        xt=[xt;S(k).X(1:end-1)'];
        yt=[yt;S(k).Y(1:end-1)'];
    end%1574701x1
    
    xt=[];
    yt=[];
    dsmooth=.1
    for k=1:length(S);
        if mod(k,100)==0,k/length(S),end
        N=length(S(k).X(1:end-1));
        c=[[0;N],[S(k).X(1:end-1);S(k).Y(1:end-1)]];
        cf=smooth_contour(c,dsmooth);
        ci=interp_contour(cf,dsmooth);
        if ~isempty(ci)
        xt=[xt;ci(1,2:end)'];
        yt=[yt;ci(2,2:end)'];
        else
          xt=[xt;mean(S(k).X(1:end-1))];
          yt=[yt;mean(S(k).Y(1:end-1))];
         
        end
    end
 
    ZPOS=XPOS+i*YPOS;
    zt=xt+i*yt;
    [NX,NY]=size(ZPOS);
   clear D
    for k=1:NX,
        if mod(k,10)==0,k/NX,end
        for j=1:NY,
            z0=xpos(j)+i*ypos(k);
            d=min(abs(z0-zt));
             D(k,j)=d;
        end
    end
         
    save -v7.3 dist_D_.mat D XPOS YPOS

    load dist_D_.mat 
    %Make High Res D by interpolating back to grid
    fl='RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc'
    lon=ncread(fl,'lon');
    lat=ncread(fl,'lat');
    z=ncread(fl,'bed_elevation');

    whos XPOS YPOS D lon lat z
    nx=length(lon);
    ny=length(lat);
    clear LON LAT
    LON=lon(:)*ones(1,ny);
    LAT=ones(nx,1)*lat(:)';
    
    whos XPOS YPOS D lon lat z LON LAT
    Di=interp2(XPOS*180/pi,YPOS*180/pi,D,LON,LAT);
   
    save -v7.3 dist_Di_.mat Di LON LAT

    

    figure
    clf;pcolor(xpos,ypos,D);shading interp;
    colorbar;axis equal;colormap('jet');hold on
    plot(real(zt),imag(zt),'k.')

    figure
    hfun=+25. - 24.5* exp(-(D/2).^2);
    hfuni=+25. - 24.5* exp(-(Di/2).^2);
    clf;pcolor(xpos,ypos,hfun);shading interp;
    colorbar;axis equal;colormap('jet');
    
    hmin=.5;%not nescesarry
    hfun(find(hfun<hmin))=hmin;
    hfuni(find(hfuni<hmin))=hmin;
        
    hfun=hfuni;      
    
    save -v7.3 hfun.mat hfun xpos ypos


    hmat.mshID = 'ELLIPSOID-GRID';
    hmat.radii = geom.radii;
    hmat.point.coord{1} = XPOS(1,:) ;
    hmat.point.coord{2} = YPOS(:,1) ;
    hmat.value = single(hfun);
    
   % opts.hfun_file='hfun'
    savemsh(opts.hfun_file,hmat) ;
    savemsh('hfunD.msh',hmat) ;

 clear hmat
    hmat.mshID = 'ELLIPSOID-GRID';
    hmat.radii = geom.radii;
    hmat.point.coord{1} = LON(1,:)*pi/180 ;
    hmat.point.coord{2} = LAT(:,1)*pi/180 ;
    hmat.value = single(hfun);
    
   % savemsh('hfunDi.msh',hmat) ;
   %for low res
    % hmat.point.coord{1} = (xpos(2:end)+xpos(1:end-1))/2;
    % hmat.point.coord{2} = (ypos(2:end)+ypos(1:end-1))/2;
    
    
     hmat.point.coord{1} = xpos;
     hmat.point.coord{2} = ypos;
     hmat.value = single(hfun);
  
   
   savemsh('hfunD.msh',hmat) ;
    %has NaNs

%------------------------------------ build mesh via JIGSAW!

    fprintf(1,'  Constructing MESH...\n');

    opts.hfun_scal = 'absolute';
    opts.hfun_hmax = +inf ;
    opts.hfun_hmin = +0.0 ;

    opts.mesh_dims = +2 ;               % 2-dim. simplexes

    opts.optm_qlim = +.95 ;
    opts.optm_iter = + 32 ;

    rbar = mean(geom.radii(:)) ;
    hbar = mean(hmat.value(:)) ;

    nlev = round(log2( ...
        rbar / sin(.4 * pi) / hbar) ) ;

    mesh = tetris(opts, nlev - 1) ;

    plotsphere(geom,mesh,hmat,topo) ;

    drawnow ;
    set(figure(1),'units','normalized', ...
    'position',[.30,.50,.25,.30]) ;

    set(figure(2),'units','normalized', ...
    'position',[.05,.50,.25,.30]) ;
    set(figure(3),'units','normalized', ...
    'position',[.05,.10,.25,.30]) ;

    savemsh('mesh_1kmUSCL_15kmOO_D',mesh) ;
    
    
mesh = loadmsh('mesh_1kmUSCL_15kmOO_D.msh');
load  GlobalCoastline5km.mat
topo=BoxSmoothTopo('RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc',3);
meshP=GlobalMesh2Planer(mesh,topo)
%meshO = CutOutLand(meshP,S,3)
meshO = CutOutLandE(meshP,S);
savemsh('mesh_1kmUSCL_15kmOO_Ocean_D.msh',meshO);


x=meshO.point.coord(:,1);
y=meshO.point.coord(:,2);
h=interp2(topo.point.coord{:,1},topo.point.coord{:,2},topo.value,x,y);
meshO.point.coord(:,3)=h(:);
e=meshO.tria3.index(:,1:3)';

close all
clf;
ph=patch_global_local(x,y,h,e,ax) 
colormap('jet');
caxis([-250,0]);
daspect([1,cos(mean(ax(3:4))*pi/180),1])
ax=[-73,-70,39,42];


'==================================================='
radii  = ones(1,3)*6371.E+00;
x=mesh.point.coord(:,1);
y=mesh.point.coord(:,2);
z=mesh.point.coord(:,3); 

S2 = R3toS2(radii,[x(:),y(:),z(:)])*180/pi;
lon=S2(:,1);
lat=S2(:,2);
h=interp2(topo.point.coord{:,1},topo.point.coord{:,2},topo.value,lon,lat);
e=mesh.tria3.index(:,1:3)';
    
ax=[-73,-70,39,42];
ph=patch_global_local(lon,lat,h,e,ax)

meshll=mesh;
meshll.point.coord(:,1)=lon(:);
meshll.point.coord(:,2)=lat(:);
meshll.point.coord(:,3)=h(:);


omesh = CutOutLandE(mesh,S)


fl='RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc'
%   [GPRJ] = project(geom,proj,'fwd');
%   [HPRJ] = project(hmat,proj,'fwd');

%mesh=loadmsh('RWPS.msh');
meshi = project(mesh,proj,'inv');

x=meshi.point.coord(:,1);y=meshi.point.coord(:,2);
z=meshi.point.coord(:,3);e=meshi.tria3.index(:,1:3)';

lont=x*180/pi;
j=find(lont>180);lon=lont;
lon(j)=lon(j)-360;lat=y*180/pi;
xe=abs(max(lon(e))-min(lon(e) ));
j=find(xe<90);
ep=e(j,:);

    
    
