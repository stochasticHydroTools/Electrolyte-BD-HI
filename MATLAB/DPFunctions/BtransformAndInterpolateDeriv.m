% Take fourier-cheb coefficients of the linear velocity back to real space,
% take its curl and interpolate to get the angular velocity of the particles.
% This uses the derivatives of the kernel to obtain angular velocity
% Input:
%   Cu,Cv,Cw - linear velocity components
%   dIx, dIy,dIz,imdIx,imdIy,imdIz - derivative interp mats and images
%   U - zeros(size(Cu))
% Output:
% uPtsDP,vPtsDP,wPtsDP - angular velocity components of the particles

function [uPtsDP,vPtsDP,wPtsDP] = BtransformAndInterpolateDeriv(Cu,Cv,Cw,dIx,dIy,dIz,...
                                                                imdIx,imdIy,imdIz,U,varargin)
    [nxy,~,Nz] = size(Cu);
    if nargin == 10
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
    % Integrate back to get the potential on the charges   
    uPtsDP = -1/2*(interpolate(dIy,wgrid,nxy,Nz)-...
                  interpolate(dIz,vgrid,nxy,Nz));
    vPtsDP = -1/2*(interpolate(dIz,ugrid,nxy,Nz)-...
                  interpolate(dIx,wgrid,nxy,Nz));
    wPtsDP = -1/2*(interpolate(dIx,vgrid,nxy,Nz)-...
                  interpolate(dIy,ugrid,nxy,Nz));
   
    if size(imdIx,1)>1 && size(imdIy,1)>1 && size(imdIz,1)>1 
      uPtsDPim = -1/2*(interpolate(imdIy,wgrid,nxy,Nz)-...
                    interpolate(imdIz,vgrid,nxy,Nz));
      vPtsDPim = -1/2*(interpolate(imdIz,ugrid,nxy,Nz)-...
                    interpolate(imdIx,wgrid,nxy,Nz));
      wPtsDPim = -1/2*(interpolate(imdIx,vgrid,nxy,Nz)-...
                    interpolate(imdIy,ugrid,nxy,Nz));
      uPtsDP = uPtsDP-uPtsDPim;
      vPtsDP = vPtsDP-vPtsDPim;
      wPtsDP = wPtsDP-wPtsDPim;
    end
end
