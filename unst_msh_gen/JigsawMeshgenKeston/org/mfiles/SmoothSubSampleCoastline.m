function [xs,ys]=SmoothSubSampleCoastline(x,y,dsmooth,lambda);
%function [xs,ys]=SmoothSubSampleCoastline(x,y,dsmooth,lambda);
%Sub sample coastline at dsmooth (m) distance. Boxcar smooth coastline to 
% dsmooth*lambda distance   
if nargin<4
    lambda=1;
end
isisland=0;

if and(x(1)==x(end),y(1)==y(end))
    isisland=1;
end
lat2m=110574.;
dx=x(2:end)-x(1:end-1);
dy=y(2:end)-y(1:end-1);
ymp=(y(2:end)+y(1:end-1))/2;
lon2m=111320.*cos(ymp*pi/180);
d=sqrt(  (dx.*lon2m).^2 + (dy.*lat2m).^2   ); 
d=[0,cumsum(d)];%distance
di=d(1):dsmooth:d(end);
N=length(di);
for k=1:N
    if mod(k,1000)==0,
        disp(['SmoothSubSampleCoastline:',num2str(k/N)]);
    end
    if isisland==1,%then use distance along a closed curve
        ind1=find(abs(d-di(k))<dsmooth*lambda);
        ind2=find(abs(d(end)-di(k)+d)<dsmooth*lambda);
        ind=union(ind1,ind2);   
    else
        ind=find(abs(d-di(k))<dsmooth*lambda);
    end
    xs(k)=mean(x(ind));
    ys(k)=mean(y(ind));
end
if isisland==1
    xs(end)=xs(1);
    ys(end)=ys(1);
end

