function omesh=RemoveNodes(mesh,j)

h.x=g.x(i);
h.y=g.y(i);
h.z=g.z(i);
A=ismember(g.e,i);
k=find(sum(A')==3)';
e=g.e(k,:);
[N,M]=size(e);

j=[];
for k=1:length(g.x);
    if mod(k,1000)==0,disp(int2str(k));end
     m=find(k==i);
     if ~isempty(m)
     j(k)=m;
 else
     j(k)=0;
 end
end
h.e=j(e);
