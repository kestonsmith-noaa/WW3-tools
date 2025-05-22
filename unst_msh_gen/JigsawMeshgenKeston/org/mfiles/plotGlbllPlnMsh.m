function plotGlbllPlnMsh(mesh,S);


x=mesh.point.coord(:,1);
y=mesh.point.coord(:,2);
z=mesh.point.coord(:,3);
e=mesh.tria3.index(:,1:3);

dx=max(x(e'))-min(x(e'));

dxmax=90;
j=find(dx<dxmax);
ep=e(j,:);
clf;pa=patch(x(ep'),y(ep'),z(ep'));shading interp;axis equal;
colormap('jet');
set(pa,'EdgeColor','k');set(pa,'EdgeAlpha',.25);
if nargin>2,
    hold on;
    for k=1:length(S),plot(S(k).X,S(k).Y,'c');end
end

ax =[  -75  -70   37 42]

axis(ax);
kprint('MeshMAET.jpg');

ax =[  -71  -70   41 42]

axis(ax);
kprint('MeshMVET.jpg');

ax =[  -71  -70.5   41.2 41.7]
axis(ax);
kprint('MeshMVhrET.jpg');
