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

    hfun = +150.0 - 112.5 * exp( -( ...
        +1.5 * (XPOS + 1.0).^2 ...
        +1.5 * (YPOS - 0.5).^2).^4) ;
    ZPOS=XPOS+i*YPOS;
    zt=xt+i*yt;
    [NX,NY]=size(ZPOS);
    hfun= 150.+0*ZPOS;

  clear D
    for k=1:NX,
        if mod(k,10)==0,k/NX,end
        for j=1:NY,
            z0=xpos(j)+i*ypos(k);
            d=min(abs(z0-zt));
             D(k,j)=d;
        end
    end
         
    save -v7.3 dist_D_.mat D

    figure
    clf;pcolor(xpos,ypos,D);shading interp;
    colorbar;axis equal;colormap('jet');hold on
    plot(real(zt),imag(zt),'k.')

    figure
    hfun=+15. - 14.* exp(-(D/2).^2);
    clf;pcolor(xpos,ypos,hfun);shading interp;
    colorbar;axis equal;colormap('jet');
    
    hmin=1;%not nescesarry
    hfun(find(hfun<hmin))=hmin;
         
    save -v7.3 hfun.mat hfun xpos ypos


    hmat.mshID = 'ELLIPSOID-GRID';
    hmat.radii = geom.radii;
    hmat.point.coord{1} = XPOS(1,:) ;
    hmat.point.coord{2} = YPOS(:,1) ;
    hmat.value = single(hfun);
    
   % opts.hfun_file='hfun'
    savemsh(opts.hfun_file,hmat) ;
   % opts.hfun_file='hfun.msh'
    

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

    savemsh('mesh_1kmUSCL_15kmOO',mesh) ;
