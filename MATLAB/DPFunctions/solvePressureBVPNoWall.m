% Solve the Poisson equation for pressure in Stokes equation
% for a doubly periodic domain with no walls (Lap(p) = div(f))

% This function only handles k>0

% Input: p_RHS - Fourier-Chebyshev coeffs of div(f).
%        BCs - returned from BCRows(Nz)
%        SIMat - second chebyshev integral matrix
%        FIMat - first chebyshev integral matrix
%        Kx,Ky - meshgrid of wave numbers
%        H - half length of z grid


% Output: Cp - Fourier-Cheb coeffs of pressure
%         Dp - Fourier-Cheb coeffs of d/dz of pressure

function [Cp,Dp,varargout] = solvePressureBVPNoWall(p_RHS,BCs,SIMat,FIMat,Kx,Ky,H)
  [nxy,~,Nz] = size(p_RHS);
  Cp = zeros(nxy,nxy,Nz);
  Dp = zeros(nxy,nxy,Nz);
  A_k = zeros(nxy,nxy,Nz,Nz);
  B_k = zeros(nxy,nxy,Nz,2);
  C_k = zeros(nxy,nxy,2,Nz);
  D_k = zeros(nxy,nxy,2,2);
  Ginv_k = zeros(nxy,nxy,2,2);
  for iX=1:nxy
    for iY=1:nxy
      kx = Kx(1,iX); ky = Ky(iY,1); k = sqrt(kx^2+ky^2);
      if k > 0
        cfhat = reshape(p_RHS(iY,iX,:),Nz,1);
        [secD, A,B,C,D,Ginv] = BVPChebInt(k,Nz,H^2*SIMat,BCs,H,cfhat,0,0); %Cp''
        Cp(iY,iX,:) = H^2*SIMat(1:Nz,:)*secD; % compute Cp
        Dp(iY,iX,:) = H*FIMat*secD; % compute derivative dCp/dz from Cp''
        A_k(iY,iX,:,:) = A;
        B_k(iY,iX,:,:) = B;
        C_k(iY,iX,:,:) = C;
        D_k(iY,iX,:,:) = D;
        Ginv_k(iY,iX,:,:) = Ginv;
      end
    end
  end
  if nargout >= 3
    varargout{1} = A_k;
  end
  if nargout >= 4
    varargout{2} = B_k;
  end
  if nargout >= 4
    varargout{3} = C_k;
  end
  if nargout >= 4
    varargout{4} = D_k;
  end
  if nargout >= 5
    varargout{5} = Ginv_k;
  end
end
