
run ~/kstart

fl='RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc'

ncdisp(fl)

lon=ncread(fl,'lon');
lat=ncread(fl,'lat');
z=ncread(fl,'bed_elevation');

whos lon lat z

imagesc(lat,lon,z')
jy=find(and(lat>25,lat<45));
jx=find(and(lon<-70,lon>-100));
zp=z(jx,jy);

clf;pcolor(lon(1:end-1),lat(1:end-1),z');shading interp

clf;pcolor(lon(jx),lat(jy),z(jx,jy)');shading interp
caxis([-6000,0])
colormap 'jet'
jz=find(and(zp'<0,zp'>-20));

lonp=lon(jx);
latp=lat(jy);
nx=length(lonp);
ny=length(latp);
LON=ones(ny,1)*lonp(:)';
LAT=latp(:)*ones(1,nx);
hold on
jz=find(and(zp'<0,zp'>-20));

plot(LON(jz),LAT(jz),'k.');

figure
clf;pcolor(LON,LAT,zp');shading interp
caxis([-6000,0])
colormap 'jet'


dx=[lon(2)-lon(1)]*111000;
dy=[lat(2)-lat(1)]*111000;
A20m=dx*dy*length(jz);
Ae=250*250/2
Ne=A20m/Ae
floor(Ne)
%east coast~ 420,000 square kilometers, 


