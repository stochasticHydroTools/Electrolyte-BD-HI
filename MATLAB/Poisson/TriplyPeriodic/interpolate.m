function pPtsDP = interpolate(I,pgrid,nxy,Nz)
    pPtsDP=I'*reshape(permute(pgrid,[2 1 3]),nxy*nxy*Nz,1);
end

