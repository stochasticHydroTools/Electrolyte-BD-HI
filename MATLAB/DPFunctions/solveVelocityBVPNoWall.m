% Solve for velocity in Stokes equation
% for a doubly periodic domain with no walls

% This function only handles k>0

% Input: u_RHS - Fourier-Cheb coeffs of (dp/dx - f_x)/eta
%        v_RHS - Fourier-Cheb coeffs of (dp/dy - f_y)/eta
%        w_RHS - Fourier-Cheb coeffs of (dp/dz - f_z)/eta
%        Cp    - Fourier-Cheb coeffs of p
%        Cf    - Fourier-Cheb coeffs of f_x
%        Cg    - Fourier-Cheb coeffs of f_y
%        BCs - returned from BCRows(Nz)
%        SIMat - second chebyshev integral matrix
%        Kx,Ky - meshgrid of wave numbers
%        Dx,Dy - Fourier 1st deriv ops
%        H - half length of z grid
%        eta - viscocity

% Output: Cu - Fourier-Cheb coeffs of x velocity
%         Cv - Fourier-Cheb coeffs of y velocity
%         Cw - Fourier-Cheb coeffs of z velocity

function [Cu,Cv,Cw] = solveVelocityBVPNoWall(u_RHS,v_RHS,w_RHS,Cp,Cf,Cg,BCs,SIMat,Kx,Ky,Dx,Dy,H,eta)
  [nxy,~,Nz] = size(u_RHS);
  % Solve the BVPs for u, v, and w
  Cu = zeros(nxy,nxy,Nz);
  Cv = zeros(nxy,nxy,Nz);
  Cw = zeros(nxy,nxy,Nz);
  for iX = 1:nxy
    for iY = 1:nxy
      cf = reshape(u_RHS(iY,iX,:),Nz,1);
      cg = reshape(v_RHS(iY,iX,:),Nz,1);
      ch = reshape(w_RHS(iY,iX,:),Nz,1);
      kx = Kx(1,iX); ky = Ky(iY,1); k = sqrt(kx^2+ky^2);
      dx = Dx(iY,iX); dy = Dy(iY,iX);
      if k > 0
        au = BVPChebInt(k,Nz,H^2*SIMat,BCs,H,cf,...
             -dx*sum(Cp(iY,iX,:))/(2*eta*k),...
             dx*sum(reshape(Cp(iY,iX,:),1,Nz).*((-1).^(0:Nz-1)))/(2*eta*k));
        cu = H^2*SIMat*au;
        
        av = BVPChebInt(k,Nz,H^2*SIMat,BCs,H,cg,...
             -dy*sum(Cp(iY,iX,:))/(2*eta*k),...
             dy*sum(reshape(Cp(iY,iX,:),1,Nz).*((-1).^(0:Nz-1)))/(2*eta*k));
        cv = H^2*SIMat*av;
        
        aw = BVPChebInt(k,Nz,H^2*SIMat,BCs,H,ch,...
             sum(Cp(iY,iX,:))/(2*eta),...
             sum(reshape(Cp(iY,iX,:),1,Nz).*((-1).^(0:Nz-1)))/(2*eta));
        cw = H^2*SIMat*aw;
        Cu(iY,iX,:) = cu;
        Cv(iY,iX,:) = cv;
        Cw(iY,iX,:) = cw;
      end
    end
  end
end
