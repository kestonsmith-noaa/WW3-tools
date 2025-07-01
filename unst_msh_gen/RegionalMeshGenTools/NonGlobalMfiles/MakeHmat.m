function MakeCoastalBoundaries1kmGlobal
% make coastlines for various smoothings of coastlines

%from example 6 aust.msh 

%geom = 
%    point: [1×1 struct]
%    edge2: [1×1 struct]
%    mshID: 'EUCLIDEAN-MESH'
%    fileV: 3

%geom.point.coord(1:10,:)
%  146.2929  -39.0150         0
%  146.2937  -39.0192         0
%  146.2846  -39.0242         0

%geom.edge2.index(1:10,:)
%           1           2           0
%           1       27577           0
%           2           3           0
           
%1km l
% 


uscl='us_coastline/tl_2023_us_coastline.shp'
US=shaperead(uscl);
NUS=length(US)
xus=[];
yus=[];
for k=1:NUS
    x=US(k).X(1:end-1);
    y=US(k).Y(1:end-1);
    [xs,ys]=SmoothSubSampleCoastlineFast(x,y,100.,25);
    xus=[xus,xs];
    yus=[yus,ys];
end

topo=BoxSmoothTopo('RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc',2);

lon = topo.point.coord{:,1};
lat = topo.point.coord{:,2};
z= topo.value;
figure;
pcolor(lon,lat,z);shading interp;
d=0*z;
nx=length(lon)
ny=length(lat)

clear d0
lat2m=110574
for j=1:ny
    j/ny
    lon2m=111320.*cos(lat(j)*pi/180);
    d0(j,:)=min (abs(  [lon(:)'-xus(:)]*lon2m + i*[ones(1,nx)*lat(j) - yus(:)]*lat2m  ) );
end
D=d0;
save -v7.3 DistToUScoast.mat lon lat D

F=exp(-D/200000);

clf;pcolor(lon,lat,F);shading interp
colormap('jet');
caxis([0,1]);
colorbar


hmat= topo.value*0+25;
 
%hmat.value = \
%        np.sqrt(np.maximum(-zlev, 0.)) / 0.5

%if z<-1000 h=25 km
z0=max(-z, 0.);
z1=max(z0,1000);
zc=-3000


z0=z;
z0(find(z>0))=0;
z0(find(z<zc))=zc;
figure;clf;pcolor(lon,lat,z0);shading interp;colorbar;colormap('jet')


H=sqrt(max(-z0, 0.));
Hc=sqrt(-zc)

Hp=25*H/Hc;
clf;pcolor(lon,lat,Hp);shading interp;colorbar;colormap('jet');

j=find(z>-20);
Hp(j)=1;
%design function
%1km at 20m 25km at 2000m
zz=0:-1:-3000;
ff=exp(-max([-20-zz]/500,0));
gg= 25-24*ff;
 clf;plot(zz,gg);

 z= topo.value;
 F=exp(-max([-20-z]/500,0));
 G= 25-24*F;
 figure;clf;pcolor(lon,lat,G);shading interp;colorbar;colormap('jet')
 %apply inside the 100km distance to US
 GP=25+0*G;
 W=exp(-D/200000);
 GPP=GP.*(1-W)+W.*G;
 figure;clf;pcolor(lon,lat,GPP);shading interp;colorbar;colormap('jet');

 hmat=topo;
 
 lon=hmat.point.coord{:,1};
 j=find(lon<90);lon(j)=180+(lon(j)+180);
 hmat.point.coord{:,1}=lon;
 hmat.value=GPP;
 savemsh('hmat_shifted',hmat);

