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

for k=1:N
    ns0(k)=length(S0(k).X(1:end-1));
end
S=S0;
save -v7.3 GlobalCoastline1km.mat S
shapewrite(S, 'GlobalCoastline1km.shp');
%this can be slow
BoundaryShape2mshVc(S,'GlobalCoastline1kmB.msh');


