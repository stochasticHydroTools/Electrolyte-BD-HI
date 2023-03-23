% Computes the spread and interpolation matrices for the
% exponential-of-semicircle (ES) kernel.
% Input:
%   xEpts, yEpts, zEpts - 1D grids in x,y and z
%   IBpts - Nx3 matrix with positions of N particles
%   iwtsz - Chebyshev weights in z
%   sup - kernel spreading width in number of points
%   phi - function handle to normalized ES kernel
% Output:
%   SWeights - spreading matrix
%   IWeights - interpolation matrix (quadrature weighted adjoint of spread)
%
% Usage:
%   The function handle phi must be instantiated as, for example:
%       w = 6; beta = 2.3*w; alpha = w*hexy/2;
%       phit = @(z) ((z/alpha).^2 <= 1).*exp(beta*(sqrt(1-(z/alpha).^2)-1));
%       normPhi = integral(@(r)phit(r),-alpha,alpha);
%       phi = @(r) phit(r)./normPhi;
function [SWeights,IWeights] = SpreadWtsES(xEpts,yEpts,zEpts,IBpts,iwtsz,sup,phi)
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
    SWeights = sparse(Nx*Ny*Nz,NIB);
    IWeights = sparse(Nx*Ny*Nz,NIB);
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
        % weights in x,y
        xwts = phi(xunwrap-IBpts(l,1));
        ywts = phi(yunwrap-IBpts(l,2));
        xywts = xwts'*ywts; xywts = xywts(:);
        % flattened indices into xy weights
        xyinds = xpts' + Nx*(ypts-1); xyinds = xyinds(:);
        % weights in z for spread and interp
        zwts = phi(zEpts(zpts)-IBpts(l,3)); 
        ziwts = zwts.*iwtsz(zpts)';
        % total weights for spread
        totwts = xywts*zwts';
        % total weights for interp
        totiwts = xywts*ziwts'*hex*hey;
        % flattened indices for spread and interp mat
        totinds = xyinds+Nx*Ny*(zpts-1);
        SWeights = SWeights+sparse(totinds(:),l,totwts(:),Nx*Ny*Nz,NIB);
        IWeights = IWeights+sparse(totinds(:),l,totiwts(:),Nx*Ny*Nz,NIB);      
    end
end
