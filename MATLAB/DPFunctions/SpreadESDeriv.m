% Computes the spread matrix for the
% derivative of the ES kernel.
% Input:
%   xEpts, yEpts, zEpts - 1D grids in x,y and z
%   IBpts - Nx3 matrix with positions of N particles
%   iwtsz - Chebyshev weights in z
%   sup - kernel spreading width in number of points
%   phi, phiDer - normalized 1D ES kernels of dipole and its derivative
%   thresh - threshold at which to truncate derivative kernel 
%   axis - chars 'x'|'y'|'z' for derivative
% Output:
%   dSWeights - deriv of ES spreading matrix
%   dIWeights - interpolation matrix (quadrature weighted adjoint of spread)

function [dSWeights,dIWeights] = SpreadESDeriv(xEpts,yEpts,zEpts,IBpts,iwtsz,sup,phi,phiDer,thresh,axis)
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
    for l=1:NIB
        % indices of zEpts within support (throwing away those outside of cutoff)
        zpts = 1:Nz;
        zpts = zpts(abs(zEpts-IBpts(l,3)) <= thresh);
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
        % throw away points outside of cutoff and correct indices
        xpts(abs(xunwrap-IBpts(l,1)) > thresh) = [];
        ypts(abs(yunwrap-IBpts(l,2)) > thresh) = [];
        xunwrap(abs(xunwrap-IBpts(l,1)) > thresh) = [];
        yunwrap(abs(yunwrap-IBpts(l,2)) > thresh) = [];
        % weights in x,y,z for partial deriv spread mat 
        if strcmp(axis,'x')
            xwts = phiDer(xunwrap-IBpts(l,1));
            ywts = phi(yunwrap-IBpts(l,2));
            zwts = phi(zEpts(zpts)-IBpts(l,3));
        elseif strcmp(axis,'y')
            xwts = phi(xunwrap-IBpts(l,1));
            ywts = phiDer(yunwrap-IBpts(l,2));
            zwts = phi(zEpts(zpts)-IBpts(l,3));
        elseif strcmp(axis,'z')
            xwts = phi(xunwrap-IBpts(l,1));
            ywts = phi(yunwrap-IBpts(l,2));
            zwts = phiDer(zEpts(zpts)-IBpts(l,3));
        else
            error('axis argument invalid - must be x,y,z');  
        end
        xywts = xwts'*ywts; xywts = xywts(:);
        % flattened indices into xy weights
        xyinds = xpts' + Nx*(ypts-1); xyinds = xyinds(:);
        % total weights for spread and interp
        totwts = xywts*zwts';
        ziwts=zwts.*iwtsz(zpts)';
        totiwts=xywts*ziwts'*hex*hey;
        % flattened indices for spread mat
        totinds = xyinds+Nx*Ny*(zpts-1);
        dSWeights = dSWeights+sparse(totinds(:),l,totwts(:),Nx*Ny*Nz,NIB);
        dIWeights = dIWeights+sparse(totinds(:),l,totiwts(:),Nx*Ny*Nz,NIB);
    end
end
