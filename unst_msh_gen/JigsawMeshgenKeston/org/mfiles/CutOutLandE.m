function omesh = CutOutLandE(mesh,S)
%load GlobalCoastline5km.mat S
plotMesh=0

if isstr(S)
    S = shaperead(S);
end

N=length(S)

x=mesh.point.coord(:,1);
y=mesh.point.coord(:,2);
z=mesh.point.coord(:,3);
e=mesh.tria3.index(:,1:3);
n1=length(x);
[ne,three]=size(e);
xe=mean(x(e'))';
ye=mean(y(e'))';
for k=1:N, ns(k)=length(S(k).X);end

jl=[];
jo=1:ne;
for k=1:N
    if mod(k,1000)==0,
        display(['progress a: ',num2str(k/N)]);,
        disp(['progress b:', num2str(  sum(ns(1:k))/sum(ns)  )]);
    end
    xl=S(k).X(1:end-1);
    yl=S(k).Y(1:end-1);
    if length(xl)>2,
        minx=min(xl);
        maxx=max(xl);
        miny=min(yl);
        maxy=max(yl);
        jx=find(and(xe>minx ,xe<maxx  ));
        jy=find(and(ye>miny ,ye<maxy  ));
        jbox=intersect(jx,jy);
        jbox=intersect(jbox,jo);% don't care about nodes out of mesh
         %ji=find(inpolygon( x(jbox),y(jbox),xl(:)',yl(:)'  ));
        ji=find(insidepoly( xe(jbox),ye(jbox),xl(:)',yl(:)'  ));
        j=jbox(ji);
        jl=[jl;j(:)];
        jo=setdiff(jo,j);
    end
end


[ne,three]=size(e);
jo=setdiff(1:ne,jl);

ep=e(jo,:);
joN=unique(ep(:));%nodes in ocean and boundary
omesh=SubsetMesh(mesh,joN);

if plotMesh==1,
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
end