function ph=patch_global_local(lon,lat,f,e,ax)

    jx=find(and(lon<ax(2),lon>ax(1)));
    jy=find(and(lat<ax(4),lat>ax(3)));
    j=intersect(jx,jy);

    A=ismember(e,j);
    kax=find(sum(A)==3)';
    eax=e(:,kax);

    ph=patch(lon(eax),lat(eax),f(eax));

