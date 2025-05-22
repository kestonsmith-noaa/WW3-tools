function geom=BoundaryShape2mshVc(S,flout)
% make coastlines for various smoothings of coastlines
%assumes all features are closed islands
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
clear geom
geom.mshID='EUCLIDEAN-MESH'
geom.fileV = 3
if isstr(S)
    S = shaperead(S);
end

N=length(S);
x0=[];
y0=[];
edge0=[];
N=length(S);
isisland=zeros(N,1);
for k=1:N
    ns(k)=length(S(k).X);
    if and(  S(k).X(end-1)==S(k).X(1) , S(k).Y(end-1)==S(k).Y(1) )
        isisland(k)=1;
    end
end

nt=sum(ns)
X=zeros(nt,1);
Y=zeros(nt,1);
n=1
for k=1:N
   l=length(S(k).X); 
   X(n:n+l-1)=S(k).X;
   Y(n:n+l-1)=S(k).Y;
   if mod(k,10000)==0,k/N,end
    n=n+l;
end

j=find(isnan(X+Y));
jg=find(~isnan(X+Y));

edge=zeros(nt-length(j),2);
edge(1:j(1)-1,1)=[1:j(1)-1]';
edge(1:j(1)-1,2)=[1:j(1)-1]'+1;
n=j(1)-1;
for k=2:N
    w=j(k)-j(k-1)-1;
    edge(n+1:n+w,1)=[j(k-1)+1:(j(k) -1)]' - k;
    edge(n+1:n+w,2)=[j(k-1)+1:(j(k) -1)]' -k + 1;
    n=j(k)-k;
end
whos X Y
 max(edge)
 min(edge)
[jj,kk]=find(edge==0)


X=X(jg);
Y=Y(jg);

close all
plot(X(edge'),Y(edge'),'k')

[ne,two]=size(edge);
nn=length(X);
geom.edge2.index=[edge,zeros(ne,1)];
geom.point.coord=[X(:),Y(:),zeros(nn,1)];
savemsh(flout,geom);
    