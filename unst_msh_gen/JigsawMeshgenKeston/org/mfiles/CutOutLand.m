function omesh = CutOutLand(mesh,S,MinWetNodes)
%load GlobalCoastline5km.mat S
if nargin < 3,MinWetNodes=3;end
if isstr(S)
    S = shaperead(S);
end

N=length(S)

x=mesh.point.coord(:,1);
y=mesh.point.coord(:,2);
z=mesh.point.coord(:,3);
e=mesh.tria3.index(:,1:3);
n1=length(x);

jl=[];
for k=1:N
    if mod(k,1000)==0,display(['progress : ',num2str(k/N)]);,end
    xl=S(k).X(1:end-1);
    yl=S(k).Y(1:end-1);
    if length(xl)>2,
        minx=min(xl);
        maxx=max(xl);
        miny=min(yl);
        maxy=max(yl);
        jx=find(and(x>minx ,x<maxx  ));
        jy=find(and(y>miny ,y<maxy  ));
        jbox=intersect(jx,jy);
        %ji=find(inpolygon( x(jbox),y(jbox),xl(:)',yl(:)'  ));
        ji=find(insidepoly( x(jbox),y(jbox),xl(:)',yl(:)'  ));
        j=jbox(ji);
        jl=[jl;j(:)];
    end
end


jo=setdiff(1:n1,jl);
[A,B]=ismember(1:n1,jo);
e0=B(e);%map to kept nodes and 0 for missing nodes
e0i=e0;
e0i(find(e0>0))=1;%index of nodes in or out

lndoci=sum(e0i');
oci=find(lndoci>=MinWetNodes);

%if IncludeEdgeE
    edgeind=find(and(lndoci>=MinWetNodes,lndoci<3));
    js=e(edgeind,:);
    js=unique(js(:));
    jo=union(jo,js);
  %  jl=setdiff(1:n1,jo);
    [A,B]=ismember(1:n1,jo);
    e0=B(e);%map to kept nodes and 0 for missing nodes
    
    e0i=e0;
    e0i(find(e0>0))=1;%index of nodes in or out

    lndoci=sum(e0i');
    oci=find(lndoci==3);
%end

  

%oci=find(lndoci>1);
clear omesh
eocean=e0(oci,:);

omesh.point.coord(:,1)=x(jo);
omesh.point.coord(:,2)=y(jo);
omesh.point.coord(:,3)=z(jo);
omesh.point.coord(:,4)=2+0*jo;
omesh.tria3.index=[eocean,zeros(length(oci),1)];


xo=omesh.point.coord(:,1);
yo=omesh.point.coord(:,2);
zo=omesh.point.coord(:,3);
eo=omesh.tria3.index(:,1:3);

figure;
ax=[-80,-65,30,45];
pa=patch(xo(eo'),yo(eo'),zo(eo'));shading interp;colorbar;caxis([-600,0]);colormap('jet');axis(ax)
set(pa,'EdgeColor','k');set(pa,'EdgeAlpha',.7);
hold on
N=length(S);
for k=1:N
    plot(S(k).X,S(k).Y,'c.-');
end