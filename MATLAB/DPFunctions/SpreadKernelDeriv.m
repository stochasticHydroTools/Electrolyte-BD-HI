% Computes the spread matrix for the
% derivative of a kernel phi. This can be used
% for both doubly and triply periodic problems

% Input:
%   xEpts, yEpts, zEpts - 1D grids in x,y and z
%   IBpts - Nx3 matrix with positions of N particles.
%           spreading is one only inside the domain.
%   iwtsz - Chebyshev weights in z
%   sup   - kernel spreading width in number of points
%   phi, phiDer - normalized 1D kernels of dipole and its derivative
%   axis - chars 'x'|'y'|'z' for spreadimg derivative d/d(axis) 
%   periodicity - of the domain (2 for doubly periodic, 3 for triple)
%               - NOTE: the code assumes Nx=Ny, hx=hy for DP
%               - NOTE: the code assumes Nx=Ny=Nz, hx=hy=hz for TP

% Output:
%   dSWeights - deriv of spreading matrix
%   dIWeights - interpolation matrix (quadrature weighted adjoint of spread)
%   See SpreadAndTransformSimple.m and BtransformAndInterpolateDeriv.m for usage details

function [dSWeights,dIWeights] = SpreadKernelDeriv(xEpts,yEpts,zEpts,IBpts,iwtsz,sup,phi,phiDer,axis,periodicity)
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
    dSWeights = sparse(Nx*Ny*Nz,NIB); dIWeights = sparse(Nx*Ny*Nz,NIB);
    ff = functions(phiDer);
    thresh = ff.workspace{1}.thresh;
    % parse axis arg 
    if strcmp(axis,'x')
        xwts_f = phiDer;
        ywts_f = phi;
        zwts_f = phi;
    elseif strcmp(axis,'y')
        xwts_f = phi;
        ywts_f = phiDer;
        zwts_f = phi;
    elseif strcmp(axis,'z')
        xwts_f = phi;
        ywts_f = phi;
        zwts_f = phiDer;
    else
        error('axis argument invalid - must be x,y,z');  
    end
    if periodicity == 2
        for l=1:NIB
            % indices of zEpts within support 
            zpts = 1:Nz;
            zpts = zpts(abs(zEpts-IBpts(l,3)) <= thresh);
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
            % unwrapped coordinates within support in x, y, z 
            xunwrap = xEpts(floorx) + mvals*hxy;
            yunwrap = yEpts(floory) + mvals*hxy;
            zunwrap = zEpts(zpts);
            % throw away points outside of cutoff and correct indices
            xpts(abs(xunwrap-IBpts(l,1)) > thresh) = []; 
            ypts(abs(yunwrap-IBpts(l,2)) > thresh) = []; 
            xunwrap(abs(xunwrap-IBpts(l,1)) > thresh) = []; 
            yunwrap(abs(yunwrap-IBpts(l,2)) > thresh) = []; 
            % weights in x,y,z for partial deriv spread mat 
            xwts = xwts_f(xunwrap-IBpts(l,1)); 
            ywts = ywts_f(yunwrap-IBpts(l,2)); 
            zwts = zwts_f(zunwrap-IBpts(l,3));
            % flattened xy weights
            xywts = xwts'*ywts; xywts = xywts(:);
            % flattened indices into xy weights
            xyinds = xpts' + Nx*(ypts-1); xyinds = xyinds(:);
            % total weights for spread and interp
            totwts = xywts*zwts';
            ziwts=zwts.*iwtsz(zpts)';
            totiwts=xywts*ziwts'*hxy*hxy;
            % flattened indices for spread mat
            totinds = xyinds+Nx*Ny*(zpts-1);
            dSWeights = dSWeights+sparse(totinds(:),l,totwts(:),Nx*Ny*Nz,NIB);
            dIWeights = dIWeights+sparse(totinds(:),l,totiwts(:),Nx*Ny*Nz,NIB);
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
            % unwrapped coordinates within support in x, y, z
            xunwrap = xEpts(floorx) + mvals*hxy;
            yunwrap = yEpts(floory) + mvals*hxy;
            zunwrap = zEpts(floorz) + mvals*hxy;
            % throw away points outside of cutoff and correct indices
            xpts(abs(xunwrap-IBpts(l,1)) > thresh) = []; 
            ypts(abs(yunwrap-IBpts(l,2)) > thresh) = []; 
            zpts(abs(zunwrap-IBpts(l,3)) > thresh) = []; 
            xunwrap(abs(xunwrap-IBpts(l,1)) > thresh) = []; 
            yunwrap(abs(yunwrap-IBpts(l,2)) > thresh) = []; 
            zpts(abs(yunwrap-IBpts(l,3)) > thresh) = []; 
            % weights in x,y,z for partial deriv spread mat 
            xwts = xwts_f(xunwrap-IBpts(l,1));
            ywts = ywts_f(yunwrap-IBpts(l,2));
            zwts = zwts_f(zunwrap-IBpts(l,3));
            xywts = xwts'*ywts; xywts = xywts(:);
            % flattened indices into xy weights
            xyinds = xpts' + Nx*(ypts-1); xyinds = xyinds(:);
            % total weights for spread and interp
            totwts = xywts*zwts';
            ziwts=zwts.*iwtsz(zpts)';
            totiwts=xywts*ziwts'*hxy*hxy;
            % flattened indices for spread mat
            totinds = xyinds+Nx*Ny*(zpts-1);
            dSWeights = dSWeights+sparse(totinds(:),l,totwts(:),Nx*Ny*Nz,NIB);
            dIWeights = dIWeights+sparse(totinds(:),l,totiwts(:),Nx*Ny*Nz,NIB);
        end
    else
      error('periodicity must be 2 or 3');
    end
end
