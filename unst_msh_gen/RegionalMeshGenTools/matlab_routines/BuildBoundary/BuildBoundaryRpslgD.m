

isplot=0;

load ../RWPSmesh/GlobalCoastline1kmUSto15km.mat

Blon=[129 20.3];
Blat=[-31,78];


earth=referenceSphere('Earth');
N=length(S);
for k=1:N
    if mod(k,1000)==0,k/N,end
    lon=S(k).X;
    j=find(lon<90);lon(j)=180+(lon(j)+180);
    S(k).X=lon;
end
lon=Blon;j=find(lon<90);lon(j)=180+(lon(j)+180);
Blon=lon;


mesh=loadmsh('RWPS.msh');
meshi = project(mesh,proj,'inv');

x=meshi.point.coord(:,1);y=meshi.point.coord(:,2);
z=meshi.point.coord(:,3);e=meshi.tria3.index(:,1:3)';

lont=x*180/pi;
j=find(lont>180);lon=lont;
lon(j)=lon(j)-360;lat=y*180/pi;
xe=abs(max(lon(e))-min(lon(e) ));
j=find(xe<90);
ep=e(j,:);



if isplot,
    
    BoundingBox(Blon,Blat,'r');
    hold on
    ax=axis;
    for k=1:N
        if mod(k,1000)==0,k/N,end
        plot(S(k).X,S(k).Y,'r');
    end
    axis(ax)
end

Bx=[Blon(1),Blon(2),Blon(2),Blon(1),Blon(1)]
By=[Blat(1),Blat(1),Blat(2),Blat(2),Blat(1)]
             
clear pslg
sxp=[];syp=[ ];
xc=[];yc=[];
N=length(S);
pslg.x=[];clf
pslg.y=[];
pslg.edges=[];
nc=0;
DXmax=720
for k=1:N
    isinbox(k)=0;
    x=S(k).X(1:end-1);% remove trailing nan (-1) and endpoint==startpoint (-2)
    y=S(k).Y(1:end-1);
        ji=find( insidepoly( x,y,Bx,By  ) );
        if length(ji)>2  
            [xi, yi,ii] = polyxpoly(x, y, Bx, By);
        %    j0=find(xi>Bx(2));xi(j0)=Bx(2);%degenerate cases
        %    j0=find(xi<Bx(1));xi(j0)=Bx(1);%degenerate cases
            xc=[xc;xi(:)];
            yc=[yc;yi(:)];
            if ~isempty(xi) %modify ob
                sxp=[sxp,x(1)];syp=[syp,y(1)];
                [iia,jja]=sort(ii(:,1));%sort to ascending order along segments
                xia=xi(jja);
                yia=yi(jja);
                iis=ii(jja,1);%sorted into ascending order along segments
                if ~insidepoly( x(1),y(1),Bx,By )% segment origonates outside box 
                    nseg=length(xia);
                    for j=1:2:nseg-1,
                        xs=[xia(j),x( iis(j)+1:iis(j+1) ),xia(j+1)];
                        ys=[yia(j),y( iis(j)+1:iis(j+1) ),yia(j+1)];
                        dx=abs(xs(2:end)-xs(1:end-1) +i*[ys(2:end)-ys(1:end-1)] );
                        if max(dx)<1,
                   
                %                j0=find(xs>Bx(2));xs(j0)=Bx(2);%degenerate cases
                %                j0=find(xs<Bx(1));xs(j0)=Bx(1);%degenerate cases

                    
                            n0=length(pslg.x);
                            pslg.x=[pslg.x,xs];
                            pslg.y=[pslg.y,ys];
                            n1=length(pslg.x);
                            nedges=[[n0+1:n1-1];[n0+2:n1]]';
                            pslg.edges=[pslg.edges;nedges];
                            nc=nc+1;pslg.chains(nc).nodes=[n0+1:n1];
                            pslg.chains(nc).index=k;
                            pslg.chains(nc).type='starts outside';
                            pslg.chains(nc).BI=1;
                        end
                    end
                else % segment starts inside outer boundary
                    nseg=length(xia);
                    xs=[xia(end),x( iis(end)+1:end ),x(1:iis(1) ),xia(1)];
                    ys=[yia(end),y( iis(end)+1:end ),y(1:iis(1) ),yia(1)];

                    dx=abs(xs(2:end)-xs(1:end-1) +i*[ys(2:end)-ys(1:end-1)] );
                    if max(dx)<1,
                    
                %                j0=find(xs>Bx(2));xs(j0)=Bx(2);%degenerate cases
                 %               j0=find(xs<Bx(1));xs(j0)=Bx(1);%degenerate cases

                        n0=length(pslg.x);
                        pslg.x=[pslg.x,xs];
                        pslg.y=[pslg.y,ys];
                        n1=length(pslg.x);
                        nedges=[[n0+1:n1-1];[n0+2:n1]]';
                        pslg.edges=[pslg.edges;nedges];
                        nc=nc+1;pslg.chains(nc).nodes=[n0+1:n1];
                        pslg.chains(nc).index=k;
                        pslg.chains(nc).type='starts inside, first seg';
                        pslg.chains(nc).BI=1;
                    end
                    for j=2:2:nseg-1,
                        xs=[xia(j),x( iis(j):iis(j+1) ),xia(j+1)];
                        ys=[yia(j),y( iis(j):iis(j+1) ),yia(j+1)];
                        
                        dx=abs(xs(2:end)-xs(1:end-1) +i*[ys(2:end)-ys(1:end-1)] );
                        if max(dx)<1,
                %                j0=find(xs>Bx(2));xs(j0)=Bx(2);%degenerate cases
                %                j0=find(xs<Bx(1));xs(j0)=Bx(1);%degenerate cases
                        
                        n0=length(pslg.x);
                        pslg.x=[pslg.x,xs];
                        pslg.y=[pslg.y,ys];
                        n1=length(pslg.x);
                        nedges=[[n0+1:n1-1];[n0+2:n1]]';
                        pslg.edges=[pslg.edges;nedges];
                        nc=nc+1;pslg.chains(nc).nodes=[n0+1:n1];
                        pslg.chains(nc).type='starts inside';
                        pslg.chains(nc).index=k;
                        pslg.chains(nc).BI=1;
                        end
                    end
                end
                
            else %isempty(xi)--> island entirely inside bounding box
                xs=[x(1:end-1)];
                ys=[y(1:end-1)];
                dx=abs(xs(2:end)-xs(1:end-1) +i*[ys(2:end)-ys(1:end-1)] );
                if and(length(xs)>2,max(dx)<1) %no degenerate islands
                    n0=length(pslg.x);
                    pslg.x=[pslg.x,xs];
                    pslg.y=[pslg.y,ys];
                    n1=length(pslg.x);
                    nedges=[[n0+1:n1-1];[n0+2:n1]]';
                    nedges=[nedges;[n1,n0+1]];
                    pslg.edges=[pslg.edges;nedges];
                    nc=nc+1;pslg.chains(nc).nodes=[n0+1:n1,n0+1];
                    pslg.chains(nc).type='interior island';
                    pslg.chains(nc).index=k;
                    pslg.chains(nc).BI=0;
                end
        end
     end
     if mod(k,1000)==0,disp(['Land segments compleate: ',num2str(k/N)]);,end
end

if isplot==1,figure;drawpslg_chainsX(pslg);title('PSLG with all box interior segments');end

%OK - now revisit outer boundary!
save -v7.3 pslgtmp.mat pslg
%Deps1=10^-8
%j=find(or( pslg.x>Bx(2)+Deps1, pslg.x<Bx(1)-Deps1) );
%pslgb=subpslg(pslg,j)
%pslg=pslgb;
%save -v7.3 pslgtmp1.mat pslg



%make outer bound
n=length(pslg.x);
pslg.x=[pslg.x,Bx(3:4)];%add northwest and north east box corner nodes
pslg.y=[pslg.y,By(3:4)];

nc=length(pslg.chains);%add northwest and north east box corner chains
pslg.chains(nc+1).nodes=n+1;
nc=length(pslg.chains);
pslg.chains(nc+1).nodes=n+2;
nc=length(pslg.chains);

xx =  153.0518
yy =  -31.0000 %south east corner of Austrilia
[m0,j0]=min(abs(pslg.x+i*pslg.y-xx-i*yy))
xx =  288.3381
yy =  -31.0000 %south west corner of S. America
[m1,j1]=min(abs(pslg.x+i*pslg.y-xx-i*yy))


nc=length(pslg.chains);
pslg.chains(nc+1).nodes=[j0,j1];
nc=length(pslg.chains);

for k=1:nc
    if mod(k,100)==0,disp(['Labeling chains compleate: ',num2str(k/nc)]);end
    spx(k)=pslg.x(pslg.chains(k).nodes(1));
    spy(k)=pslg.y(pslg.chains(k).nodes(1));
    epx(k)=pslg.x(pslg.chains(k).nodes(end));
    epy(k)=pslg.y(pslg.chains(k).nodes(end));
 end
 
%make boundary order index along boundary for start points
SN=10.*[max(abs(pslg.x))+max(abs(pslg.y))];%large number to seperate edges
c=0*spy;
Deps1=10^-10
j=find(abs(spy-By(1))<Deps1);
c(j)=SN+spx(j);
j=find(abs(spx-Bx(2))<Deps1);
c(j)=2*SN+spy(j);
j=find(abs(spy-By(3))<Deps1);
c(j)=3*SN-spx(j);
j=find(abs(spx-Bx(4))<Deps1);
c(j)=4*SN-spy(j);

%make boundary order index along boundary for end points
d=0*spy;
j=find(abs(epy-By(1))<Deps1);
d(j)=1*SN+epx(j);
j=find(abs(epx-Bx(2))<Deps1);
d(j)=2*SN+epy(j);
j=find(abs(epy-By(3))<Deps1);
d(j)=3*SN-epx(j);
j=find(abs(epx-Bx(4))<Deps1);
d(j)=4*SN-epy(j);

figure;plot(pslg.x,pslg.y,'k.');hold on
epi=j0;%start at lower left boundary corner
obc=[nc];%start at lower left boundary corner
v=SN+pslg.x(epi);
j=find(c>v);
[mm,m]=min(c(j));
l=j(m);
epi=[epi,pslg.chains(l).nodes];
obc=[obc,l];
n=0;
while(epi(end)~=epi(1))
    v=d(l);
    j=find(c>v);
    [mm,m]=min(c(j));
    l=j(m);
    epi=[epi,pslg.chains(l).nodes];
    obc=[obc,l];
    plot(pslg.x(epi),pslg.y(epi),'r.-');pause(.000001)
        
end

% add edges in obc
for k=1:length(epi)-1
    pslg.edges=[pslg.edges;[epi(k),epi(k+1)]];
end

xob=pslg.x(epi);
yob=pslg.y(epi);
nodelist=[];
for k=1:nc
    if mod(k,100)==0,disp(['Finding interior chains compleate: ',num2str(k/nc)]);,end
    if pslg.chains(k).nodes(1)==pslg.chains(k).nodes(end)
        n=pslg.chains(k).nodes(1);
        if inside( pslg.x(n), pslg.y(n),xob,yob)>0,
            nodelist=union(nodelist,pslg.chains(k).nodes);
        end
    end
end
nodelist=union(nodelist,epi);
nn=length(pslg.x);
pslgb=subpslg(pslg,setdiff(1:nn,nodelist));

if(0)
    %remove chains in Outer Boundary Chain (OBC) from interior chains
    ic=setdiff(1:nc,obc);
    pslg.chains=pslg.chains(ic);

    %add outer boundary chain to chain list
    nc=length(pslg.chains)
    pslg.chains(nc+1).nodes=epi;
    nc=length(pslg.chains);

    %remove chains and nodes outside OBC
    xob=pslg.x(epi);
    yob=pslg.y(epi);
    %[inpoly,onpoly]=insidepoly( pslg.x,pslg.y,xob,yob,'tol',1./10000 ); %misses some points
    %[inpoly,onpoly]=insidepoly( pslg.x,pslg.y,xob,yob );
    [inpoly]=inside( pslg.x,pslg.y,xob,yob );
    j=find(inpoly>0);
    nn=length(pslg.x);
    pslgb=subpslg(pslg,setdiff(1:nn,j));
end

if isplot==1
    figure;drawpslg(pslgb);title('PSLG with all points interior to final outer boundary')
end

%remove "random" duplicate nodes
z=pslgb.x+i*pslgb.y;
[zu,j,k]=unique(z);
nn=length(pslgb.x);
pslgb=subpslg(pslgb,setdiff(1:nn,j));

%remove duplicate edges
pslgb.edgesS=sort(pslgb.edges')';
pslgb.edgesSU=unique(pslgb.edgesS,'rows')
pslgb.edges=pslgb.edgesSU;

pslg=pslgb

save PSLGboundary.mat pslg

geom=pslg2geom(pslg)

savemsh('PSLGboundary.msh',geom)
