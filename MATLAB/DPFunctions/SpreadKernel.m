% Computes the spread and interpolation matrices given a 
% normalized spreading kernel phi. This can be used
% for both doubly and triply periodic problems

% Input:
%   xEpts, yEpts, zEpts - 1D grids in x,y and z
%   IBpts - Nx3 matrix with positions of N particles.
%           spreading is one only inside the domain.
%   iwtsz - Chebyshev weights in z
%   sup - kernel spreading width in number of points
%   phi - function handle to normalized kernel
%   periodicity - of the domain (2 for doubly periodic, 3 for triple)
%               - NOTE: the code assumes Nx=Ny, hx=hy for DP
%               - NOTE: the code assumes Nx=Ny=Nz, hx=hy=hz for TP
% Output:
%   SWeights - spreading matrix
%              the force on grid is given by SWeights * forces
%   IWeights - interpolation matrix (quadrature weighted adjoint of spread)
%              particle velocity is given by IWeights * Ugrid,
%              where Ugrid is fluid velocity on the grid
%   See SpreadAndTransformSimple.m and BtransformAndInterpolate.m for more details


function [SWeights,IWeights] = SpreadKernel(xEpts,yEpts,zEpts,IBpts,iwtsz,sup,phi,periodicity)
    Nx = length(xEpts);
    Ny = length(yEpts);
    Nz = length(zEpts);
    [NIB,~] = size(IBpts);
    hxy = xEpts(2)-xEpts(1);
    aex = min(xEpts);
    aey = min(yEpts);
    aez = min(zEpts);
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
    if periodicity == 2
        for l=1:NIB
            % indices of zEpts within support
            zpts = 1:Nz;
            zpts = zpts(abs(zEpts-IBpts(l,3)) <= sup/2*hxy);
            % index of grid pt to left of particle center in x
            floorx = mod(floor((IBpts(l,1)-aex)/hxy),Nx) + 1;
            % correct mvals for when particle is closer to right face of cell i
            if oddsup && abs(xEpts(floorx) - IBpts(l,1)) > hxy/2
                floorx = floorx + 1;
            end      
            % periodic indices within support in x
            xpts = floorx-down:floorx+up;
            xpts = mod(xpts,Nx);
            xpts(xpts==0) = Nx;
            % index of grid pt to left of particle center in y
            floory = mod(floor((IBpts(l,2)-aey)/hxy),Ny) + 1;
            % correct mvals for when particle is closer to right face of cell i
            if oddsup && abs(yEpts(floory) - IBpts(l,2)) > hxy/2
                floory = floory + 1;
            end  
            % periodic indices within support in y
            ypts = floory-down:floory+up;
            ypts = mod(ypts,Ny);
            ypts(ypts==0) = Ny;
            % unwrapped coordinates within support in x and y
            xunwrap = xEpts(floorx) + mvals*hxy;
            yunwrap = yEpts(floory) + mvals*hxy;
            % weights in x,y
            xwts = phi(xunwrap-IBpts(l,1));
            ywts = phi(yunwrap-IBpts(l,2));
            xywts = xwts'*ywts; xywts = xywts(:);            
            zwts = phi(zEpts(zpts)-IBpts(l,3)); 
            % flattened indices into xy weights
            xyinds = xpts' + Nx*(ypts-1); xyinds = xyinds(:);
            % weights in z for spread and interp
            ziwts = zwts.*iwtsz(zpts)';
            % total weights for spread
            totwts = xywts*zwts';
            % total weights for interp
            totiwts = xywts*ziwts'*hxy*hxy;
            % flattened indices for spread and interp mat
            totinds = xyinds+Nx*Ny*(zpts-1);
            SWeights = SWeights+sparse(totinds(:),l,totwts(:),Nx*Ny*Nz,NIB);
            IWeights = IWeights+sparse(totinds(:),l,totiwts(:),Nx*Ny*Nz,NIB);      
        end
    elseif periodicity == 3
        for l=1:NIB
            floorz = mod(floor((IBpts(l,3)-aez)/hxy),Nz) + 1;
            % correct mvals for when particle is closer to right face of cell i
            if oddsup && abs(zEpts(floorz) - IBpts(l,3)) > hxy/2
                floorz = floorz + 1;
            end         
            % periodic indices within support in z
            zpts = floorz-down:floorz+up;
            zpts = mod(zpts,Nz);
            zpts(zpts==0) = Nz;
            
            % index of grid pt to left of particle center in x
            floorx = mod(floor((IBpts(l,1)-aex)/hxy),Nx) + 1;
            % correct mvals for when particle is closer to right face of cell i
            if oddsup && abs(xEpts(floorx) - IBpts(l,1)) > hxy/2
                floorx = floorx + 1;
            end      
            % periodic indices within support in x
            xpts = floorx-down:floorx+up;
            xpts = mod(xpts,Nx);
            xpts(xpts==0) = Nx;
            % index of grid pt to left of particle center in y
            floory = mod(floor((IBpts(l,2)-aey)/hxy),Ny) + 1;
            % correct mvals for when particle is closer to right face of cell i
            if oddsup && abs(yEpts(floory) - IBpts(l,2)) > hxy/2
                floory = floory + 1;
            end  
            % periodic indices within support in y
            ypts = floory-down:floory+up;
            ypts = mod(ypts,Ny);
            ypts(ypts==0) = Ny;
            % unwrapped coordinates within support in x and y
            xunwrap = xEpts(floorx) + mvals*hxy;
            yunwrap = yEpts(floory) + mvals*hxy;
            zunwrap = zEpts(floorz) + mvals*hxy;
            % weights in x,y
            xwts = phi(xunwrap-IBpts(l,1));
            ywts = phi(yunwrap-IBpts(l,2));
            xywts = xwts'*ywts; xywts = xywts(:);
            % flattened indices into xy weights
            xyinds = xpts' + Nx*(ypts-1); xyinds = xyinds(:);
            % weights in z for spread and interp
            zwts = phi(zunwrap'-IBpts(l,3)); 
            ziwts = zwts * hxy;
            % total weights for spread
            totwts = xywts*zwts';
            % total weights for interp
            totiwts = xywts*ziwts'*hxy*hxy;
            % flattened indices for spread and interp mat
            totinds = xyinds+Nx*Ny*(zpts-1);
            SWeights = SWeights+sparse(totinds(:),l,totwts(:),Nx*Ny*Nz,NIB);
            IWeights = IWeights+sparse(totinds(:),l,totiwts(:),Nx*Ny*Nz,NIB);      
        end
    else
      error('periodicity must be 2 or 3');
    end
end
