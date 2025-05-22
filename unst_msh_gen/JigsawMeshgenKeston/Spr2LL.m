M=loadmsh('RWPS.Sphr.msh');

p=loadmsh('PSLGboundary5km.msh')
xmin=min(p.point.coord(:,1))
xmax=max(p.point.coord(:,1))
ymin=min(p.point.coord(:,2))
ymax=max(p.point.coord(:,2))


proj.prjID = 'stereographic'
proj.radii = +6.371E+003
proj.xbase = +0.500 * (xmin + xmax) * pi / 180.
proj.ybase = +0.500 * (ymin + ymax) * pi / 180.

Mp = project(M,proj,'inv')
Mp.point.coord(:,1:2)=Mp.point.coord(:,1:2)*180/pi;

figure;
plot(Mp.point.coord(:,1),Mp.point.coord(:,2),'k.',p.point.coord(:,1),p.point.coord(:,2),'r.')
savemsh('RWPS.LL.msh',Mp);

topo=loadmsh('topo_shifted.msh');

x=Mp.point.coord(:,1);
y=Mp.point.coord(:,2);
xt=double(topo.point.coord{:,1});
yt=double(topo.point.coord{:,2});

z=interp2(xt,yt, double(topo.value),x,y);

xp=x-360;
A=[xp(:),y(:),z(:)];
Mp.point.coord(:,1:3)=A;
%Mp.point.coord(:,3)=z(:);

Mp.point=rmfield(Mp.point,'power')
savemsh('RWPS.LL.H.msh',Mp);
e=Mp.tria3.index(:,1:3);

figure;ph=patch(x(e'),y(e'),z(e'));shading interp
set(ph,'EdgeColor','k');
set(ph,'EdgeAlpha',.2);

ax=[-73,-70,39,42];
ax(1:2)=ax(1:2)+360;
axis(ax)

hold on;plot(p.point.coord(:,1),p.point.coord(:,2),'b.');