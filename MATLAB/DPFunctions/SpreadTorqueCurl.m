% This function spreads the curl of torques on the grid
% while taking into account the images
 
% Input: 
%   dSx,dSy,dSz - spread matrices for kernel derivatives
%   torques   - torques on particles (Np x 3)
%   imdSx,imdSy,imdSz  - image spread matrices for kernel derivatives
%   imTorques - image torques
%             - these should be EQUAL to the torque on the corresponding particle
%             - we end up subtracting the spread of the image torque curl from that
%               of the torque curl on the actual particles
%   nxy,Nz   - number of points in x/y and z

% Output:
%   curlx,curly,curlz - x,y and z component of spread curl of torques on grid

function [curlx,curly,curlz,dxT,dyT,dzT] = SpreadTorqueCurl(dSx,dSy,dSz,torques,imdSx,imdSy,imdSz,imTorques,nxy,Nz)
    % spread curl of torques
    [curlx,curly,curlz,dxT,dyT,dzT] = spreadCurl(dSx,dSy,dSz,torques,nxy,Nz);
    [nIm,~] = size(imTorques);
    % apply near wall correction on torque curl if needed
    if nIm ~= 0
        [imCurlx,imCurly,imCurlz,imdxT, imdyT, imdzT] = spreadCurl(imdSx,imdSy,imdSz,imTorques,nxy,Nz);
        curlx = curlx-imCurlx; curly = curly-imCurly; curlz = curlz-imCurlz;
        dxT = dxT - imdxT; dyT = dyT - imdyT; dzT = dzT - imdzT;
    end
end
