function [gridfx, gridfy, gridfz] = spread(S,forces,nxy,Nz)
    [~,dim]=size(forces);
    gridf=full(S*forces);
    gridfx=permute(reshape(gridf(:,1),nxy,nxy,Nz),[2 1 3]);
    if (dim == 3)
        gridfy=permute(reshape(gridf(:,2),nxy,nxy,Nz),[2 1 3]);
        gridfz=permute(reshape(gridf(:,3),nxy,nxy,Nz),[2 1 3]);
    end
end