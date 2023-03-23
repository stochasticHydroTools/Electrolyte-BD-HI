% Take fourier-cheb coefficients  of linear or angular velocity back to real space
% and interpolate to get the lin/ang velocity of the particles.
% This does not use the derivatives of the kernel to obtain angular velocity
% Input:
%   Cu,Cv,Cw - either lin or ang velocity components
%   I, imI - interp mat and the image
%   U - zeros(size(Cu))
% Output:
% uPtsDP,vPtsDP,wPtsDP - linear or angular velocity components of the particles

function [uPtsDP,vPtsDP,wPtsDP] = BtransformAndInterpolate(Cu,Cv,Cw,I,imI,U,varargin)
    [nxy,~,Nz] = size(Cu);
    if nargin == 6
      ugrid = btransform(Cu,Nz,U); 
      vgrid = btransform(Cv,Nz,U); 
      wgrid = btransform(Cw,Nz,U);
    else
      Tu = varargin{1};
      ugrid = btransform(Cu,Nz,U,Tu);
      vgrid = btransform(Cv,Nz,U,Tu);
      wgrid = btransform(Cw,Nz,U,Tu);
      [Nz,~] = size(Tu); 
    end
    % interpolate velocity on particles
    uPtsDP = interpolate(I,ugrid,nxy,Nz); 
    vPtsDP = interpolate(I,vgrid,nxy,Nz);
    wPtsDP = interpolate(I,wgrid,nxy,Nz);
    % interpolate velocity of images and correct
    if size(imI,1)>1
      uPtsDPim = interpolate(imI,ugrid,nxy,Nz); 
      vPtsDPim = interpolate(imI,vgrid,nxy,Nz);
      wPtsDPim = interpolate(imI,wgrid,nxy,Nz);
      uPtsDP = uPtsDP-uPtsDPim; 
      vPtsDP = vPtsDP-vPtsDPim; 
      wPtsDP = wPtsDP-wPtsDPim;
    end 
end
