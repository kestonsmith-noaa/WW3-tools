

flcl='/mnt/sdb/keston/RWPSmesh/AliMeshTools/ne_110m_coastline.shp'

 load /mnt/sdb/keston/RWPSmesh/AliMeshTools/Global10km.mat 


 flcl='/mnt/sdb/keston/RWPSmesh/AliMeshTools/coastlines-split-4326/lines.shp'
 [A,R] = shaperead(flcl);

 ns=length(A);
%hold on
%for k=1:ns,
%    if mod(k,1000)==0,k/ns,end
%    plot(A(k).X,A(k).Y,'k');
%end
 lonB= [129.117, 15.7099]
 latB=[ -30.1389, 76.4989]

 latEUS=[18,45.556]
 lonEUS=[-63.973,-98.1185]

 jx=find(and( lon>min(lonEUS),lon<max(lonEUS) ));
 jy=find(and( lat>min(latEUS),lat<max(latEUS) ));
zp=z(jx,jy);
jz=find(and(zp<0,zp>-20));


 N=length(A);
p.x=[]
p.y=[]
p.edges=[]
p.holes.x=[]
p.holes.y=[]
clear qi;
j=0;
for k=1:N
    a=A(k);
    ynm=nanmean(a.Y);
    xnm=nanmean(a.X);
    if mod(k,1000)==0,k/N,dx,end   
        
    if(and(  ynm>0 , and(  xnm>min(lon),xnm<max(lon)  ))),
 %       if xnm>-50,dx=10000;end %europe /africa
 %       if xnm<-50,dx=1000;end %N. America
        if max(a.X)>-50,dx=10000/100000;end %europe /africa
        if max(a.X)<=-50,dx=1000/100000;end %N. America
        q.x=a.X(1:end-1);
        q.y=a.Y(1:end-1);
        if q.x(end)==q.x(1),%island
            holeflag=1;
            q.x=q.x(1:end-1);
            q.y=q.y(1:end-1);
        else
             holeflag=0;
        end
        n=length(q.x);
        q.edges=[[1:n-1]',[2:n]'];
        if holeflag==1,
            q.edges=[q.edges;[n,1]];
        end
        c=pslg2contour(q);
        cs=smooth_contour(c,2*dx);
        ci=interp_contour(cs,dx);
        j=j+1;
        qi(j)=contour2pslg(ci,holeflag);
        
        %p=joinpslg_v2(qi,p);
    end
end
save -v7.3 NorthAtlPSLG.mat qi


initjig;                            % load jigsaw
