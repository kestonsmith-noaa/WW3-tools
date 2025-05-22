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
           
%1km 
clear 
geom.mshID='EUCLIDEAN-MESH'
geom.fileV = 3

earth=referenceSphere('Earth')

gcfl='./GlobalCoast/GSHHS_shp/f/GSHHS_f_L1.shp'
S = shaperead(gcfl);
N=length(S);
isisland=zeros(N,1);
for k=1:N
    ns(k)=length(S(k).X(1:end-1));
    if and(  S(k).X(end-1)==S(k).X(1) , S(k).Y(end-1)==S(k).Y(1) )
        isisland(k)=1;
    end
end
if sum(isisland)==N
    disp(['All features in ',gcfl,' are closed islands'])
else
    disp([int2str(sum(isisland)),'  features in ',gcfl,...
        ' are closed islands. out of:',int2str(N),' total features'])
end

S0=S;
n=0;
for k=1:N
    x=S(k).X(1:end-1);
    y=S(k).Y(1:end-1);
    [xs,ys]=SmoothSubSampleCoastlineFast(x,y,50.,10);
    if length(xs)>2,
        n=n+1;
        S0(n).X=[xs(:);NaN]';
        S0(n).Y=[ys(:);NaN]';
        S0(n).area = areaint(ys,xs,earth) / 10^6;
    end
            
    if mod(k,1000)==0,
        disp(['progress a:', num2str(  sum(ns(1:k))/sum(ns)  )]);
        disp(['progress b:', num2str(  k/N  )]);
    end
end

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

S0=S;

load GlobalCoastline15km.mat %for speed
N=length(S);
for k=1:N
    ns(k)=length(S(k).X(1:end-1));
end
 S=rmfield(S,'dUS')
earth=referenceSphere('Earth')
for k=1:N
    x=S(k).X(1:end-1);
    y=S(k).Y(1:end-1);

  %  [xs,ys]=SmoothSubSampleCoastlineFast(x,y,100.,25);
    for j=1:length(x);
        S(k).dUS(j)=min(distance(y(j),x(j),yus,xus,earth));
    end
    if mod(k,100)==0,
        disp(['progress a:', num2str(  sum(ns(1:k))/sum(ns)  )]);
        disp(['progress b:', num2str(  k/N  )]);
    end
end

clear dUSmin
for k=1:N
    m=min(S(k).dUS);
    if ~isempty(m)
        dUSmin(k)=m;
    else
        dUSmin(k)=NaN;
    end
end

clear dUSmean
for k=1:N
    m=mean(S(k).dUS);
    if ~isempty(m)
        dUSmean(k)=m;
    else
        dUSmean(k)=NaN;
    end
end


save -v7.3 GlobalCoastline15km.mat S dUSmin dUSmean%with distance to US coast

contUSlon=[ -125.0 , -66.9 ]
contUSlat=[ 24.5  , 49.0]
Alaskalon=[	-180.,-130.]
Alaskalat=[	51.229088,71.352561]

cutoff_f= .06;

close all
M=colormap(jet(64));
for k=1:N

    x=S(k).X(1:end-1);
    y=S(k).Y(1:end-1);
    if ~isnan(dUSmin(k))
        rgb=[dUSmean(k)-min(dUSmean)]./[max(dUSmean)-min(dUSmean)];
        c=interp1([0:63]'/63,M,rgb);
        plot(x,y, 'color', c);hold on
    end
    if mod(k,100)==0,
        disp(['progress a:', num2str(  sum(ns(1:k))/sum(ns)  )]);
        disp(['progress b:', num2str(  k/N  )]);
    end
end

f=[dUSmean-min(dUSmean)]./[max(dUSmean)-min(dUSmean)];

cutoff_f= .06;
figure;
plot(xus,yus, 'k.');hold on
M=colormap(jet(64));
for k=1:N

    x=S(k).X(1:end-1);
    y=S(k).Y(1:end-1);
    if f(k)<.1
        rgb=f(k);
        c=interp1([0:63]'/63,M,rgb);
        plot(x,y, 'color', c);
    end
    if mod(k,100)==0,
        disp(['progress a:', num2str(  sum(ns(1:k))/sum(ns)  )]);
        disp(['progress b:', num2str(  k/N  )]);
    end
end
colormap(M)
colorbar

 save -v7.3 GlobalCoastline15km.mat S dUSmin dUSmean f cutoff_f


 close all
 clear
load GlobalCoastline15km.mat
S15=S;
load GlobalCoastline1km.mat
S1=S;
load GlobalCoastline5km.mat
S5=S;
gcfl='./GlobalCoast/GSHHS_shp/f/GSHHS_f_L1.shp'
S = shaperead(gcfl)

id=S(1).id
source=S(1).source
Geometry=S(1).Geometry
level=S(1).level
parent_id=S(1).parent_id
sibling_id= S(1).sibling_id

N=length(S1)
figure;hist(f,100)
usthresh=.06;
nearthresh=.1;
n=0;
earth=referenceSphere('Earth');
clear S
for k=1:N
    if f(k)<usthresh,
        x=S1(k).X(1:end-1);
        y=S1(k).Y(1:end-1);
        if length(x)>2,
            n=n+1;
            S(n).area = areaint(y,x,earth) / 10^6;
            S(n).X=S1(k).X;
            S(n).Y=S1(k).Y;
            S(n).BoundingBox=[[min(x);max(x)],[min(y);max(y)]];

            S(n).id=id;
            S(n).source=source;
            S(n).Geometry=Geometry;
            S(n).level=level;
            S(n).parent_id=parent_id;            
            S(n).sibling_id=sibling_id;

        end

    elseif f(k)<nearthresh
        x=S5(k).X(1:end-1);
        y=S5(k).Y(1:end-1);
        if length(x)>2,
            n=n+1;
            S(n).area = areaint(y,x,earth) / 10^6;
            S(n).X=S5(k).X;
            S(n).Y=S5(k).Y;
            S(n).BoundingBox=[[min(x);max(x)],[min(y);max(y)]];
             S(n).id=id;
            S(n).source=source;
            S(n).Geometry=Geometry;
            S(n).level=level;
            S(n).parent_id=parent_id;            
            S(n).sibling_id=sibling_id;
      end
    else
        
        x=S15(k).X(1:end-1);
        y=S15(k).Y(1:end-1);
        if length(x)>2,
            n=n+1;
            S(n).area = areaint(y,x,earth) / 10^6;
            S(n).X=S15(k).X;
            S(n).Y=S15(k).Y;
            S(n).BoundingBox=[[min(x);max(x)],[min(y);max(y)]];
            S(n).id=id;
            S(n).source=source;
            S(n).Geometry=Geometry;
            S(n).level=level;
            S(n).parent_id=parent_id;            
            S(n).sibling_id=sibling_id;
       end
    end
    if mod(k,1000)==0,
        disp(['progress b:', num2str(  k/N  )]);
    end
end



save -v7.3 GlobalCoastline1kmUSto15km.mat S
shapewrite(S, 'GlobalCoastline1kmUSto15km.shp');
%this can be slow
BoundaryShape2mshVc(S,'GlobalCoastline1kmUSto15km.msh');


