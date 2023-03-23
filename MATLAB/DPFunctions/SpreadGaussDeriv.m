% Computes the spread matrix for the
% derivative of a Gaussian kernel.
% Input:
%   xEpts, yEpts, zEpts - 1D grids in x,y and z
%   IBpts - Nx3 matrix with positions of N particles
%   iwtsz - Chebyshev weights in z
%   sup - kernel spreading width in number of points
%   gw - width of Gaussian
%   axis - chars 'x'|'y'|'z'|'zy'|'zx'|'zz' for derivative
% Output:
%   dSWeights - deriv of Gaussian spreading matrix
%   IWeights - interpolation matrix (quadrature weighted adjoint of spread)

function [dSWeights, dIWeights] = SpreadGausDeriv(xEpts,yEpts,zEpts,IBpts,iwtsz,sup,gw,axis)
    Nx = length(xEpts);
    Ny = length(yEpts);
    Nz = length(zEpts);
    [NIB,~] = size(IBpts);
    hex = xEpts(2)-xEpts(1);
    hey = yEpts(2)-yEpts(1);
    aex = min(xEpts);
    aey = min(yEpts);
    down = sup/2-1;
    up = sup/2;
    oddsup = false;
    if (mod(sup,2) == 1)
        down = floor(sup/2);
        up = down;
        oddsup = true;
    end
    mvals = -down:up;
    dSWeights = sparse(Nx*Ny*Nz,NIB); dIWeights = sparse(Nx*Ny*Nz,NIB);
    % normalization for Gaussian
    normFact = 1/sqrt(8*pi^3*gw^6);
    for l=1:NIB
        % indices of zEpts within support
        zpts = 1:Nz;
        zpts = zpts(abs(zEpts-IBpts(l,3)) <= sup/2*hex);
        % index of grid pt to left of particle center in x
        floorx = mod(floor((IBpts(l,1)-aex)/hex),Nx) + 1;
        % correct mvals for when particle is closer to right face of cell i
        if oddsup && abs(xEpts(floorx) - IBpts(l,1)) > hex/2
            floorx = floorx + 1;
        end      
        % periodic indices within support in x
        xpts = floorx-down:floorx+up;
        xpts = mod(xpts,Nx);
        xpts(xpts==0) = Nx;
        % index of grid pt to left of particle center in y
        floory = mod(floor((IBpts(l,2)-aey)/hey),Ny) + 1;
        % correct mvals for when particle is closer to right face of cell i
        if oddsup && abs(yEpts(floory) - IBpts(l,2)) > hey/2
            floory = floory + 1;
        end  
        % periodic indices within support in y
        ypts = floory-down:floory+up;
        ypts = mod(ypts,Ny);
        ypts(ypts==0) = Ny;
        % unwrapped coordinates within support in x and y
        xunwrap = xEpts(floorx) + mvals*hex;
        yunwrap = yEpts(floory) + mvals*hey;
        % weights in x,y,z for spread
        if strcmp(axis,'x')
            xwts = -(xunwrap-IBpts(l,1))./gw^2.*exp(-(xunwrap-IBpts(l,1)).^2./(2*gw^2));
            ywts = exp(-(yunwrap-IBpts(l,2)).^2/(2*gw^2));
            zwts = exp(-(zEpts(zpts)-IBpts(l,3)).^2/(2*gw^2));
        elseif strcmp(axis,'y')
            xwts = exp(-(xunwrap-IBpts(l,1)).^2/(2*gw^2));
            ywts = -(yunwrap-IBpts(l,2))./gw^2.*exp(-(yunwrap-IBpts(l,2)).^2./(2*gw^2));
            zwts = exp(-(zEpts(zpts)-IBpts(l,3)).^2/(2*gw^2));
        elseif strcmp(axis,'z')
            xwts = exp(-(xunwrap-IBpts(l,1)).^2/(2*gw^2));
            ywts = exp(-(yunwrap-IBpts(l,2)).^2/(2*gw^2));
            zwts = -(zEpts(zpts)-IBpts(l,3))./gw^2.*exp(-(zEpts(zpts)-IBpts(l,3)).^2./(2*gw^2)); 
        elseif strcmp(axis,'zx')
            xwts = (xunwrap-IBpts(l,1)).*exp(-(xunwrap-IBpts(l,1)).^2./(2*gw^2));
            ywts = exp(-(yunwrap-IBpts(l,2)).^2/(2*gw^2));
            zwts = (zEpts(zpts)-IBpts(l,3))./gw^4.*exp(-(zEpts(zpts)-IBpts(l,3)).^2./(2*gw^2));
        
        elseif strcmp(axis,'zy')
            xwts = exp(-(xunwrap-IBpts(l,1)).^2./(2*gw^2));
            ywts = (yunwrap-IBpts(l,2)).*exp(-(yunwrap-IBpts(l,2)).^2/(2*gw^2));
            zwts = (zEpts(zpts)-IBpts(l,3))./gw^4.*exp(-(zEpts(zpts)-IBpts(l,3)).^2./(2*gw^2));
        elseif strcmp(axis,'zz')
            xwts = exp(-(xunwrap-IBpts(l,1)).^2./(2*gw^2));
            ywts = exp(-(yunwrap-IBpts(l,2)).^2/(2*gw^2));
            zwts = exp(-(zEpts(zpts)-IBpts(l,3)).^2./(2*gw^2)).*((zEpts(zpts)-IBpts(l,3)).^2./gw^4 - 1/gw^2);
        else
            error('axis argument invalid - must be x,y,z,zx,zy or zz');  
        end
        xywts = xwts'*ywts; xywts = xywts(:);
        % flattened indices into xy weights
        xyinds = xpts' + Nx*(ypts-1); xyinds = xyinds(:);
        % total weights for spread and interp
        totwts = normFact*xywts*zwts';
        ziwts=zwts.*iwtsz(zpts)';
        totiwts=normFact*xywts*ziwts'*hex*hey;
        % flattened indices for spread mat
        totinds = xyinds+Nx*Ny*(zpts-1);
        dSWeights = dSWeights+sparse(totinds(:),l,totwts(:),Nx*Ny*Nz,NIB);
        dIWeights = dIWeights+sparse(totinds(:),l,totiwts(:),Nx*Ny*Nz,NIB);
    end
end
