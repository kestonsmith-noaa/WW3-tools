function MakeCoastalBoundaries
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
for k=1:N
    ns(k)=length(S(k).X(1:end-1));
end
for k=1:N
    x=S(k).X(1:end-1);
    y=S(k).Y(1:end-1);

    [xs,ys]=SmoothSubSampleCoastlineFast(x,y,50.,5);
    S0(k).X=[xs(:);NaN]';
    S0(k).Y=[ys(:);NaN]';
    if mod(k,10)==0,
        disp(['progress:', num2str(  sum(ns(1:k))/sum(ns)  )]);
    end
end
N0=length(S0);
for k=1:N
    ns0(k)=length(S0(k).X(1:end-1));
end

S=S0;
save -v7.3 GlobalCoastline500m.mat S
BoundaryShape2msh(S,'GlobalCoastline500m.msh')

