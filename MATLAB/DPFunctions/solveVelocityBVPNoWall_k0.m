% Solve for velocity in Stokes equation
% for a doubly periodic domain with no walls

% This function only handles k=0

% Input: u_RHS - Fourier-Cheb coeffs of (dp/dx - f_x)/eta for k=0
%        v_RHS - Fourier-Cheb coeffs of (dp/dy - f_y)/eta for k=0
%        w_RHS - Fourier-Cheb coeffs of (dp/dz - f_z)/eta for k=0
%        BCs - returned from BCRows(Nz)
%        SIMat - second chebyshev integral matrix
%        H - half length of z grid
%        Cf - Fourier-Cheb coeffs of f_x
%        Cg - Fourier-Cheb coeffs f f_y
%        uvints - precomputed integrals of cheb poly for vel correction
%        eta - viscocity
%        k0 - switch to handle the k=0 mode the k=0 mode 
%           - if k0=0, the k=0 mode of the RHS for velocity will be 0
%             and the solutions are also 0
%           - if k0=1, the k=0 mode of the RHS for velocity will not be 0,
%             the solutions will not be 0, and they will be corrected

% Output: Cu - Fourier-Cheb coeffs of x velocity for k=0
%         Cv - Fourier-Cheb coeffs of y velocity for k=0
%         Cw - Fourier-Cheb coeffs of z velocity for k=0

function [Cu,Cv,Cw] = solveVelocityBVPNoWall_k0(u_RHS,v_RHS,w_RHS,BCs,SIMat,H,Cf,Cg,uvints,eta,k0)
  % sol will be 0 if k0 = 0
  if k0 == 0
    if sum(u_RHS(1,1,:)) ~= 0 || sum(v_RHS(1,1,:)) ~= 0 || sum(w_RHS(1,1,:)) ~= 0
      error('velocity RHS for k=0 is not 0')
    end 
    Cu = 0; Cv = 0; Cw = 0;
  % if k=0 of RHS is not 0, we solve with homogenous BCs and
  % add the linear in z term for the null space   
  elseif k0 == 1
    [nxy,~,Nz] = size(u_RHS);
    cf = reshape(u_RHS(1,1,:),Nz,1);
    cg = reshape(v_RHS(1,1,:),Nz,1);
    ch = reshape(w_RHS(1,1,:),Nz,1);
    au = BVPChebInt(0,Nz,H^2*SIMat,BCs,H,cf,0,0);
    av = BVPChebInt(0,Nz,H^2*SIMat,BCs,H,cg,0,0);
    aw = BVPChebInt(0,Nz,H^2*SIMat,BCs,H,ch,0,0);
    Cu = H^2*SIMat*au;
    Cv = H^2*SIMat*av;
    Cw = H^2*SIMat*aw;
    % Add the correction to the first mode of u and v. For w, Az=0. (see eqs 29,39 of Ondrej's report)
    Cu(2) = Cu(2)+0.5/eta*uvints*reshape(Cf(1,1,:),Nz,1);
    Cv(2) = Cv(2)+0.5/eta*uvints*reshape(Cg(1,1,:),Nz,1);
  end
end
