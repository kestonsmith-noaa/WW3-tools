function MakeCoastalBoundaries5kmGlobal
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
 
gcfl='./GlobalCoast/GSHHS_shp/f/GSHHS_f_L1.shp'
S = shaperead(gcfl);
N=length(S);
S0=S;
for k=1:N
    ns(k)=length(S(k).X(1:end-1));
end
for k=1:N
    x=S(k).X(1:end-1);
    y=S(k).Y(1:end-1);
    [xs,ys]=SmoothSubSampleCoastlineFast(x,y,100.,25);
    S0(k).X=[xs(:);NaN]';
    S0(k).Y=[ys(:);NaN]';
    if mod(k,100)==0,
        disp(['progress a:', num2str(  sum(ns(1:k))/sum(ns)  )]);
        disp(['progress b:', num2str(  k/N  )]);
    end
end
for k=1:N
    ns0(k)=length(S0(k).X(1:end-1));
end
S=S0;
save -v7.3 GlobalCoastline5km.mat S
BoundaryShape2msh(S,'GlobalCoastline5km.msh');

gcfl='./GlobalCoast/GSHHS_shp/f/GSHHS_f_L1.shp'
Sorig = shaperead(gcfl);
load GlobalCoastline5km.mat %S
N=length(S);
for k=1:N
    if mod(k,1000)==0,k/N,end
    S(k).Geometry=Sorig(k).Geometry;
    S(k).BoundingBox=[[min(S(k).X);max(S(k).X)],[min(S(k).Y);max(S(k).Y)]];
    S(k).id=Sorig(k).id;
    S(k).level = Sorig(k).level; 
    S(k).source = Sorig(k).source; 
    S(k).parent_id = Sorig(k).parent_id ;
    S(k).sibling_id = Sorig(k).sibling_id; 
    S(k).area = Sorig(k).area;
end
shapewrite(S, 'GlobalCoastline5km.shp');



