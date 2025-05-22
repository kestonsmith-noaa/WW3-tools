
fl='/mnt/sdb/keston/GlobalBathymetry/GEBCO_2024_sub_ice_topo.nc'
z=ncread(fl,'elevation');
w=ones(3,3)/9;
z=conv2(z,w,'same');
z=z(1:3:end,1:3:end);
save -v7.3 GEBCOss3.mat z


lon=ncread(fl,'lon');
lat=ncread(fl,'lat');
z=ncread(fl,'elevation');
w=ones(3,3)/9;
z=conv2(z,w,'same');
z=z(1:3:end,1:3:end);
lon=lon(1:3:end);
lat=lat(1:3:end);
save -v7.3 GEBCOss3.mat lon lat z

clear


jx=find(and(lon<15,lon>-102));
jy=find(lat>=0);
whos z lon lat jx jy

z=z(jx,jy);
lon=lon(jx);
lat=lat(jy);

%close all;pcolor(lon,lat,z');shading interp;axis equal
save -v7.3 NorthAtl.mat lon lat z

W3=ones(3,3)/9;
z3=conv2(z,W3,'same');
lon3=lon(1:3:end);
lat3=lat(1:3:end);
z3=z3(1:3:end,1:3:end);

lon=lon3;lat=lat3;z=z3;
save -v7.3 NorthAtlSS3.mat lon lat z

close all;pcolor(lon,lat,z');shading interp;axis equal;colormap('jet');
caxis([-7000,0]);

flcl='/mnt/sda/keston/GlobalCoastline/GlobalCoastStreet/coastlines-split-4326/lines.shp'
[A,R] = shaperead(flcl);

figure
N=length(A);
for k=1:N
    if mod(k,1000)==0,k/N,end
    a=A(k);
    ynm=nanmean(a.Y);
    if ynm>0
        xnm=nanmean(a.X);
        if(and(xnm>min(lon),xnm<max(lon)))
            plot(a.X,a.Y,'k');hold on;
        end
    end
end

%coast line does not have lakes- ! need right file for ths

close all
flw='/mnt/sda/keston/GlobalCoastline/simplified-water-polygons-split-3857/simplified_water_polygons.shp'
[B,Rb] = shaperead(flw);
N=length(B);
clear SL
for k=1:N
    if mod(k,1000)==0,k/N,end
    a=B(k);
    SL(k)=length(a.Y);
    ynm=nanmean(a.Y);
 %   if ynm>0
        xnm=nanmean(a.X);
%        if(and(xnm>min(lon),xnm<max(lon)))
            plot(a.X,a.Y,'k');hold on
%        end
%    end
end

close all
N=length(Rb);
for k=1:N
    if mod(k,1000)==0,k/N,end
    a=Rb(k);
    ynm=nanmean(a.x);
    if ynm>0
        xnm=nanmean(a.x);
        if(and(xnm>min(lon),xnm<max(lon)))
            plot(a.x,a.y,'k');hold on
        end
    end
end


%Make Pslg Europe at 10km US at 1KM
clear
close all

flcl='/mnt/sda/keston/GlobalCoastline/GlobalCoastStreet/coastlines-split-4326/lines.shp'
[A,R] = shaperead(flcl);


%Length in km of 1° of latitude = always 111.32 km
%Length in km of 1° of longitude = 40075 km * cos( latitude ) / 360

figure

N=length(A);
p.x=[]
p.y=[]
p.edges=[]
p.holes.x=[]
p.holes.y=[]
clear qi;
j=0;
for k=1:N
    a=A(k);
    ynm=nanmean(a.Y);
    xnm=nanmean(a.X);
    if mod(k,1000)==0,k/N,dx,end   
        
    if(and(  ynm>0 , and(  xnm>min(lon),xnm<max(lon)  ))),
 %       if xnm>-50,dx=10000;end %europe /africa
 %       if xnm<-50,dx=1000;end %N. America
        if max(a.X)>-50,dx=10000/100000;end %europe /africa
        if max(a.X)<=-50,dx=1000/100000;end %N. America
        q.x=a.X(1:end-1);
        q.y=a.Y(1:end-1);
        if q.x(end)==q.x(1),%island
            holeflag=1;
            q.x=q.x(1:end-1);
            q.y=q.y(1:end-1);
        else
             holeflag=0;
        end
        n=length(q.x);
        q.edges=[[1:n-1]',[2:n]'];
        if holeflag==1,
            q.edges=[q.edges;[n,1]];
        end
        c=pslg2contour(q);
        cs=smooth_contour(c,2*dx);
        ci=interp_contour(cs,dx);
        j=j+1;
        qi(j)=contour2pslg(ci,holeflag);
        
        %p=joinpslg_v2(qi,p);
    end
end
save -v7.3 NorthAtlPSLG.mat qi

        
nf=length(qi)
for k=1:nf
     if mod(k,1000)==0,k/nf,dx,end   
   ns(k)=length(qi(k).x);
    z=qi(k).x+i*qi(k).y;
    d(k)=sum(abs(z(2:end)-z(1:end-1)));
end

figure;
k=1;
plot(qi(k).x,qi(k).y,'k.');hold on
for k=2:nf
   
    if mod(k,1000)==0,k/nf,end   
    %e=qi(k).edges';
    %plot(qi(k).x(e),qi(k).y(e),'k.-');
    plot(qi(k).x,qi(k).y,'k.-');
end

p.x=[]
p.y=[]
p.edges=[]
p.holes.x=[]
p.holes.y=[]

for k=1:nf
    if mod(k,1000)==0,k/nf,end   
    p=joinpslg_v2(qi(k),p);
end
save -v7.3 NorthAtlPSLGLL.mat qi p

figure
drawpslg(p)


N=length(A);
for k=1:N
    bs(k)=length(A(k).X);
end


N=length(A);
p.x=[]
p.y=[]
p.edges=[]
p.holes.x=[]
p.holes.y=[]
clear qi;
j=0;
for k=1:N
    a=A(k);
    ynm=nanmean(a.Y);
    xnm=nanmean(a.X);
    if mod(k,1000)==0,k/N,dx,end   
        
    if(and(  ynm>0 , and(  xnm>min(lon),xnm<max(lon)  ))),
 %       if xnm>-50,dx=10000;end %europe /africa
 %       if xnm<-50,dx=1000;end %N. America
        if max(a.X)>-50,dx=20000;end %europe /africa
        if max(a.X)<=-50,dx=2000;end %N. America
        q.x=a.X(1:end-1);
        q.y=a.Y(1:end-1);
        if q.x(end)==q.x(1),%island
            holeflag=1;
            q.x=q.x(1:end-1);
            q.y=q.y(1:end-1);
        else
             holeflag=0;
        end
        n=length(q.x);
        q.edges=[[1:n-1]',[2:n]'];
        if holeflag==1,
            q.edges=[q.edges;[n,1]];
        end
        c=pslg2contour(q);
        cs=smooth_contourLLu(c,2*dx);
        ci=interp_contourLLu(cs,dx);
        j=j+1;
        qi(j)=contour2pslg(ci,holeflag);
        
        %p=joinpslg_v2(qi,p);
    end
end


save -v7.3 NorthAtlPSLGLLu.mat qi




     clear ns d  
nf=length(qi)
for k=1:nf
     if mod(k,1000)==0,k/nf,dx,end   
   ns(k)=length(qi(k).x);
    z=qi(k).x+i*qi(k).y;
    d(k)=sum(abs(z(2:end)-z(1:end-1)));
end

p.x=[]
p.y=[]
p.edges=[]
p.holes.x=[]
p.holes.y=[]

for k=1:nf
    if mod(k,1000)==0,k/nf,end   
   if ns(k)>2,p=joinpslg_v2(qi(k),p);end
end
save -v7.3 NorthAtlPSLGLLu.mat qi p


