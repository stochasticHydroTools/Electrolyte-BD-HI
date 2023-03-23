% Solve the Poisson equation for pressure in Stokes equation
% for a doubly periodic domain with no walls (Lap(p) = div(f))

% This function only handles k=0

% Input: p_RHS - Fourier-Chebyshev coeffs of div(f).
%        BCs - returned from BCRows(Nz)
%        SIMat - second chebyshev integral matrix
%        FIMat - first chebyshev integral matrix
%        H - half length of z grid
%        Ch    - Fourier-Cheb coeffs of f_z
%        pints - precomputed integrals of cheb poly for pressure correction
%        k0 - switch to handle the k=0 mode the k=0 mode 
%           - if k0=0, the k=0 mode of the RHS for pressure will be 0
%             and the solutions are also 0
%           - if k0=1, the k=0 mode of the RHS for pressure will not be 0,
%             the solutions will not be 0, and they will be corrected

% Output: Cp - Fourier-Cheb coeffs of pressure for k=0
%         Dp - Fourier-Cheb coeffs of d/dz of pressure for k=0

function [Cp,Dp] = solvePressureBVPNoWall_k0(p_RHS,BCs,SIMat,FIMat,H,Ch,pints,k0)
  % sol will be 0 if k0 = 0
  if k0 == 0
    if p_RHS(1,1,:) ~= 0
      error('pressure RHS for k=0 is not 0');
    end
    Cp = 0; Dp = 0; 
  % if k=0 of RHS is not 0, we solve the pressure poblem with homogenous BCs
  % and must correct the 2nd Chebyshev coefficient of that mode (see eq. 18 in Ondrej's report)
  elseif k0 == 1
    [nxy,~,Nz] = size(p_RHS);
    cfhat = reshape(p_RHS(1,1,:),Nz,1);
    secD = BVPChebInt(0,Nz,H^2*SIMat,BCs,H,cfhat,0,0); %Cp''
    Cp = H^2*SIMat(1:Nz,:)*secD; % compute Cp
    Dp = H*FIMat*secD; % compute derivative dCp/dz from Cp''
    % Correction to pressure solution (add linear in z term for null space)
    Cp(2) = Cp(2)+0.5*pints*reshape(Ch(1,1,:),Nz,1);
  end 
end
