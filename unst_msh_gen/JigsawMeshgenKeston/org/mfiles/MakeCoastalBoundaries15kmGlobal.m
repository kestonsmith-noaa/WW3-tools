function MakeCoastalBoundaries15kmGlobal
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
S0=S;
for k=1:N
    ns(k)=length(S(k).X(1:end-1));
end
for k=1:N
    x=S(k).X(1:end-1);
    y=S(k).Y(1:end-1);
    [xs,ys]=SmoothSubSampleCoastlineFast(x,y,100.,75);
    S0(k).X=[xs(:);NaN]';
    S0(k).Y=[ys(:);NaN]';
    if mod(k,1000)==0,
        disp(['progress a:', num2str(  sum(ns(1:k))/sum(ns)  )]);
        disp(['progress b:', num2str(  k/N  )]);
    end
    S0(k).area = areaint(S(k).Y(1:end-1),S(k).X(1:end-1),earth) / 10^6;
end
for k=1:N
    ns0(k)=length(S0(k).X(1:end-1));
end
S=S0;
save -v7.3 GlobalCoastline15km.mat S
BoundaryShape2msh(S,'GlobalCoastline15km.msh');
shapewrite(S, 'GlobalCoastline15km.shp');



