% This function computes the spread and interp matrices
% as well as the ones for images for spreading curl of torques onto the grid
% That is, the matrices are for the derivatives of the dipole kernel

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

%   phid - normalized 1D kernel handle for dipole
%        - One can use the following:
%               1) phid = ES(alpham,betam);
%               2) phid = Gaussian(sigma,alpham);

%   phid_der - normalized 1D derivative kernel handle for dipole
%            - One can use the following:
%               1) phid_der = ES_d(alpham,betam);
%               2) phid_der = Gaussian_d(sigma,alpham);

% Output:
%   dSx,dSy,dSz,imdSx,imdSy,imdSz - derivative spreading matrix and the one for images
%   dIx,dIy,dIz,imdIx,imdIy,imdIz - derivative interpolation matrix and the one for images

function [dSx,dIx,dSy,dIy,dSz,dIz,imdSx,imdIx,imdSy,imdIy,imdSz,imdIz] = ...
    SpreadInterpMatRot(xE,yE,zpts,zwts,wd,phid,phid_der,pts,imPtst)
    
    % derivs of dipole
    [dSx,dIx] = SpreadKernelDeriv(xE,yE,zpts,pts,zwts,wd,phid,phid_der,'x',2);
    [dSy,dIy] = SpreadKernelDeriv(xE,yE,zpts,pts,zwts,wd,phid,phid_der,'y',2);
    [dSz,dIz] = SpreadKernelDeriv(xE,yE,zpts,pts,zwts,wd,phid,phid_der,'z',2);
    [nIm,~] = size(imPtst);
    if nIm ~= 0
        [imdSx,imdIx] = SpreadKernelDeriv(xE,yE,zpts,imPtst,zwts,wd,phid,phid_der,'x',2);
        [imdSy,imdIy] = SpreadKernelDeriv(xE,yE,zpts,imPtst,zwts,wd,phid,phid_der,'y',2);
        [imdSz,imdIz] = SpreadKernelDeriv(xE,yE,zpts,imPtst,zwts,wd,phid,phid_der,'z',2);
    else
        imdSx = 0; imdIx = 0;
        imdSy = 0; imdIy = 0;
        imdSz = 0; imdIz = 0;
    end
end
