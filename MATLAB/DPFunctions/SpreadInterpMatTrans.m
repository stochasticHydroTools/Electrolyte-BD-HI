% This function computes the spread and interp matrices
% as well as the ones for images for spreading forces onto the grid

% Input:
%   xE, yE, zpts - 1D grids in x,y and z 
%                - usually, x,y \in [-L,L-h] for h the xy grid spacing, z \in [a,b].
%                - That is, the x,y grids should not include L as an endpoint 
%   zwts - clenshaw curtis weights

%   pts - Nx3 matrix with positions of N particles
%       - spreading is done only inside the domain.
%       - any part of the kernel out of the domain is handled via images

%   imPts - positions of image particles

%   wm - spreading width for monopole kernel

%   phim - normalized 1D kernel handle for monopole
%        - One can use the following:
%               1) phim = ES(alpham,betam);
%               2) phim = Gaussian(sigma,alpham);
%               3) phim = @(r) 1/h*stnd4pt(r/h);
%               4) phim = @(r) 1/h*flex5pt(r/h,(38-sqrt(69))/60);
%               5) phim = @(r) 1/h*flex6pt(r/h,59/60-sqrt(29)/20);

% Output:
%   S,imS - spreading matrix and the one for images
%   I,imI - interpolation matrix and the one for images

function [S,I,imS,imI] = SpreadInterpMatTrans(xE,yE,zpts,zwts,wm,phim,pts,imPts)
    [S,I] = SpreadKernel(xE,yE,zpts,pts,zwts,wm,phim,2);
    [nIm,~] = size(imPts);
    if nIm ~= 0
        [imS,imI] = SpreadKernel(xE,yE,zpts,imPts,zwts,wm,phim,2); % force image
    else
        imS = 0; imI = 0;
    end 
end
