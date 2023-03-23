function gridfx = spread(S,forces,N)
    gridf=full(S*forces);
    gridfx=permute(reshape(gridf(:,1),nxy,nxy,Nz),[2 1 3]);
end
