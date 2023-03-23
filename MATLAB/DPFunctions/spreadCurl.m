% spreads the curl of the torque using the kernel's
% analytical 1st derivatives
function [curlx, curly, curlz, dxT, dyT, dzT] = spreadCurl(dSx,dSy,dSz,torques,nxy,Nz)
    dxT = full(dSx*torques);
    dyT = full(dSy*torques);
    dzT = full(dSz*torques);
    %disp([sum(abs(dxT(:))), sum(abs(dyT(:))), sum(abs(dzT(:)))])
    dyTz = permute(reshape(dyT(:,3),nxy,nxy,Nz),[2 1 3]);
    dzTy = permute(reshape(dzT(:,2),nxy,nxy,Nz),[2 1 3]);
    dzTx = permute(reshape(dzT(:,1),nxy,nxy,Nz),[2 1 3]);
    dxTz = permute(reshape(dxT(:,3),nxy,nxy,Nz),[2 1 3]);
    dxTy = permute(reshape(dxT(:,2),nxy,nxy,Nz),[2 1 3]);
    dyTx = permute(reshape(dyT(:,1),nxy,nxy,Nz),[2 1 3]);
    curlx = dyTz - dzTy;
    curly = dzTx - dxTz;
    curlz = dxTy - dyTx;
end
