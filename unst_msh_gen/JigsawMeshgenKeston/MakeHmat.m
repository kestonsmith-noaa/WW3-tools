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
clear
close all

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
    DLON= mod(lon(:)'-xus(:),360) ;% large matrix
%    d0(j,:)=min (abs(  [lon(:)'-xus(:)]*lon2m + i*[ones(1,nx)*lat(j) - yus(:)]*lat2m  ) );
    ld1=min (  abs(  [ DLON ]*lon2m + i*[ones(1,nx)*lat(j) - yus(:)]*lat2m  ) );
    ld2=min (  abs(  [360 - DLON ]*lon2m + i*[ones(1,nx)*lat(j) - yus(:)]*lat2m  ) ); 360-ld1;
    d0(j,:)=min(  [ ld1(:),ld2(:) ]'  );
end
D=d0;
save -v7.3 DistToUScoast.mat lon lat D

 lambda=111000*4.;
 W=exp(-D/lambda);
 %W2=exp(-[D/lambda].^2);
 
 close all
figure;clf;pcolor(lon,lat,W);shading interp
colormap('jet');
caxis([0,1]);
colorbar

% G2= 15-14*W2;
% figure;clf;pcolor(lon,lat,G2);shading interp;colorbar;colormap('jet');
% title('G2 exp(-[D/lambda].^2)');
 

close all
BGR=25.
SHR=.5
G= BGR-(BGR-SHR)*W;
 
 hmat=topo;
 
 lon=hmat.point.coord{:,1};
 j=find(lon<90);lon(j)=180+(lon(j)+180);
 j0=setdiff(1:length(lon),j);j0=j0(:);
 lon=lon([j0(:);j(:)]);
 G1=[G(:,j0),G(:,j)];
 hmat.point.coord{:,1}=lon;
 hmat.value=G1;

 figure;clf;pcolor(hmat.point.coord{:,1},hmat.point.coord{:,2},hmat.value);
shading interp;colorbar;colormap('jet')
 


 savemsh('hfun_shifted',hmat);

%Now interpolate to resolution of topography
clear topo
topo=topo2msh('RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc');
lon=topo.point.coord{:,1};
j=find(lon<90);lon(j)=180+(lon(j)+180);
topo.point.coord{:,1}=lon;
% fix break at IDL
jn=find(lon<180);
jp=find(lon>=180);
lon=lon(:);
lon=[lon(jn);lon(jp)];
topo.point.coord{:,1}=lon;
topo.value=[topo.value(jn,:);topo.value(jp,:)]';

savemsh('topo_shifted',topo);

x=topo.point.coord{:,1};
y=topo.point.coord{:,2};
nx=length(x);
ny=length(y);
X=ones(ny,1)*[x(:)'];
Y=y(:)*ones(1,nx);
Gi=interp2(hmat.point.coord{:,1},hmat.point.coord{:,2},hmat.value,X,Y);
topo.value=Gi;
topo.value(find( isnan(topo.value ) ) )=max(hmat.value(:));
savemsh('hfun_shiftedHR',topo);
