% Computes the components to assemble the RHS of transformed Stokes equations in 
% DOUBLY PERIODIC domains with 0,1 or 2 walls. We also return the interpolation 
% matrices to recover particle velocity (lin/rot) using a spreading kernel.

% Input:
%   forces - on the particles

%   xE, yE, zpts - 1D grids in x,y and z 
%                - usually, x,y \in [-L,L-h] for h the xy grid spacing, z \in [a,b].
%                - That is, the x,y grids should not include L as an endpoint 
%   zwts - clenshaw curtis weights

%   pts - Nx3 matrix with positions of N particles
%       - spreading is done only inside the domain.
%       - any part of the kernel out of the domain is handled via images

%   wm - spreading width for monopole kernel

%   phiMnorm - normalized 1D kernel handle for monopole
%            - One can use the following:
%               1) phiMnorm = ES(alpham,betam);
%               2) phiMnorm = Gaussian(sigma,alpham);
%               3) phiMnorm = @(r) 1/h*stnd4pt(r/h);
%               4) phiMnorm = @(r) 1/h*flex5pt(r/h,(38-sqrt(69))/60);
%               5) phiMnorm = @(r) 1/h*flex6pt(r/h,59/60-sqrt(29)/20);

%   domType - 0 for no walls
%           - 1 for bottom wall
%           - 2 for bottom and top wall

%   msep = [msep_down,msep_up] - maximum distance between down/up z boundaries and particle
%          for monopole to get an image (or exit with error). See setup_grid.m:
%        - If domType = 0 and a particle is within msep_down of lower bndry 
%          or msep_up of upper bndry, the function will exit with an error. 
%          That is, if there are no walls, no particles should be within msep 
%          of either upper or lower boundaries.
%          Suggested setting is msep_down = msep_up = 1.5*wh/2
%        - If domType = 1 and particle iP is such that (pts(iP,3)-a) < msep(1),
%          an image about z = a will be used for its monopole. 
%          If (b-pts(iP,3)) < msep(2), the function will exit with an error.
%          That is, no particles should be within msep(2) of the upper bndry.
%          Suggested setting is msep_down = wh/2, msep_up - 1.5*wh/2
%        - If domType = 2 and particle iP is such that (pts(iP,3)-a) < msep(1),
%          an image about z = a will be used for its monopole.
%          Similarly, if (b-pts(iP,3)) < msep(2), an image
%          about z = b will be used.
%          Suggested setting is msep_down = msep_up = wh/2

%   varargin{1} = torques - on the particles

%   varargin{2} = wd - spreading width for dipole kernel

%   varargin{3} = phiDnorm - normalized 1D kernel handle for dipole
%                          - See phiMnorm doc above for options

%   varargin{4} = phiDnorm_d - handle to derivative of 1D kernel for dipole
%                            - One can use the following
%                              1) phiDnorm_d = ES_d(alphad,betad);
%                              2) phiDnorm_d = Gaussian_d(sigma,alpham);

%   varargin{5} = dsep - maximum distance between boundaries and particle
%                        for dipole to get an image
%                      - Note, msep >= dsep

% Output :
%   Cf,Cg,Ch,Df,Dg,Dh - Components of RHS of stokes equations
%                     - Each of these are Fourier-Chebyshev coefficients
%                       of a spread quantity on the grid.
%                     - The code below computes (but does not return)
%                       spreading matrices S and imS for all particles, 
%                       where imS handles the images. 
%                     - For example, the code spreads forces
%                       as (S * forces - imS * forces), so that forces on the z-grid
%                       boundary are 0.
%                     - Cf,Cg,Ch are the x,y,z components of the RHS of Stokes 
%                       equations given in eq. 1.6 of Sachin's report.
%                     - Df,Dg,Dh are derivatives in z of Cf,Cg and Ch
%                     - Eg) For domType = 0, the RHS to the pressure Poisson equation
%                       is assembled as p_RHS = Dx.*Cf + Dy.*Cg + Dh, where Dx and 
%                       Dy are meshgrids of Fourier 1st derivative ops (1i*Kx,1i*Ky).
%                       The RHS for the velocity BVP is given by  
%                       u_RHS = (Dx.*Cp-Cf)/eta; v_RHS = (Dy.*Cp-Cg)/eta; w_RHS = (Dp-Ch)/eta;,
%                       where Cp and Dp are the coefficients of the pressure solution and its derivative.
%                       See eqs 4 and 17-19 of Ondrej's report 
%                     - Eg) For domType = 1|2, the RHS to auxiliary Stokes equation
%                       is assembled for each k as RHS = [sparse(zeros(Nz,1));df-dx*ch;dg-dy*ch],
%                       where df = Df(iY,iX,:), dg = Dg(iY,iX,:), ch = Ch(iY,iX,:), 
%                       dx = Dx(1,iX), dy = Dy(iY,1). See kSolveBVP[One|Two]Wall.m
%                       and eqs. 1.22-1.24 of Sachin's report

%   I,imI, - monopole interp matrices for all particles (I) and also for
%            images about the wall(s) for particles whose 
%            kernel overlaps the wall (imI). imI = 0 if domType = 0.
%          - To obtain the an interpolated quantity for the particles
%            use, for example: Upts = (I * Ugrid - imI * Ugrid). This is done
%            in BtransformAndInterpolate.m

%   varargout{1} = dIx - dipole x-derivative interp matrix
%   varargout{2} = dIy - dipole y-derivative interp matrix
%   varargout{3} = dIz - dipole z-derivative interp matrix
%   varargout{4} = imdIx - image interp matrix for x-deriv of dipole (=0 if domType = 0)
%   varargout{5} = imdIy - image interp matrix for y-deriv of dipole (=0 if domType = 0)
%   varargout{6} = imdIz - image interp matrix for z-deriv of dipole (=0 if domType = 0)

% Usage:

% For translation-Only problems, call
%   For domType = 0 (note, imI will be 0):
%     [Cf,Cg,Ch,~,~,Dh,I,imI] = ...
%       SpreadAndTransformSimple(forces,xE,yE,zpts,zwts,...
%                                pts,wm,phiMnorm,0,msep);

%   For domType = 1|2:
%     [Cf,Cg,Ch,Df,Dg,~,I,imI] = ...
%       SpreadAndTransformSimple(forces,xE,yE,zpts,zwts,...
%                                pts,wm,phiMnorm,1|2,msep);

%   After solving for Fourier-Chebyshev coefficients of velocity Cu,Cv,Cw,
%   interpolate particle velocity while taking into account the images with
%   [uPts,vPts,wPts] = BtransformAndInterpolate(Cu,Cv,Cw,I,imI,U);

% For problems with rotation, call
%   For no wall domType = 0: (note, imI,imdIx,imdIy,imdIz will be 0)
%     [Cf,Cg,Ch,~,~,Dh,I,imI,dIx,dIy,dIz,imdIx,imdIy,imdIz] = ...
%       SpreadAndTransformSimple(forces,xE,yE,zpts,zwts,...
%                                pts,wm,phiMnorm,0,msep,...
%                                torques,wd,phiDnorm,phiDnorm_d,dsep);

%   For one or two wall domType = 1|2:
%     [Cf,Cg,Ch,Df,Dg,~,I,imI,dIx,dIy,dIz,imdIx,imdIy,imdIz] = ...
%       SpreadAndTransformSimple(forces,xE,yE,zpts,zwts,...
%                                pts,wm,phiMnorm,1|2,msep,...
%                                torques,wd,phiDnorm,phiDnorm_d,dsep);

%   After solving for Fourier-Chebyshev coefficients of velocity Cu,Cv,Cw, 
%   take curl of linear fluid velocity while interpolating and get the
%   angular velocity of the particles with
%   [uPtst,vPtst,wPtst] = BtransformAndInterpolateDeriv(Cu,Cv,Cw,dIx,dIy,dIz,...
%                                                       imdIx,imdIy,imdIz,U);


% NOTE: For translation only, there should be 10 input and 8 output args
%       For trans-rot, there should be 15 input and at least 14 output args


function [Cf,Cg,Ch,Df,Dg,Dh,I,imI,varargout] = ...
    SpreadAndTransformSimple(forces,xE,yE,zpts,zwts,pts,...
                             wm,phiMnorm,domType,msep,varargin)
    
    % enforce validity of input list and parse
    if nargin > 10
        if (nargin < 15)
            error('not enough inputs for torque');
        else
            torques = varargin{1};
            wd = varargin{2};
            phiDnorm = varargin{3};
            phiDnorm_d = varargin{4};
            dsep = varargin{5};
            has_torque = true;
        end
    else
        has_torque = false;
    end

    % get grid info
    nxy = length(xE);
    Nz = length(zpts);
    b = zpts(1); a = zpts(end);  
    H = (b-a)/2;
    [Np,~] = size(pts);

    % define particle images and forces/torques based on msep/dsep
    % also check validity of particle configuration given domType
    iIm = 0; iImt = 0; 
    imPts = []; imPtst = []; 
    imForces = []; imTorques = [];
    for iP = 1:Np
        if (pts(iP,3)-a) < msep(1) && abs((pts(iP,3)-a)-msep(1)) > 1e-15
            if domType == 0
                error(strcat('Particle iP=', num2str(iP),' is too close to z=',num2str(a),' for DP no wall'));
            end
            iIm = iIm+1;
            imPts(iIm,:) = pts(iP,:); imPts(iIm,3) = 2*a - pts(iP,3);
            imForces(iIm,:) = 1.0*forces(iP,:);
            if has_torque && (pts(iP,3)-a) < dsep(1) && abs((pts(iP,3)-a)-dsep(1)) > 1e-15
                iImt = iImt+1;
                imPtst(iImt,:) = imPts(iIm,:);
                imTorques(iImt,:) = 1.0*torques(iP,:);
            end
        elseif (b-pts(iP,3)) < msep(2) && abs((b-pts(iP,3))-msep(2)) > 1e-15
            if domType == 0
                error(strcat('Particle iP=', num2str(iP),' is too close to z=',num2str(b),' for DP no wall'));
            elseif domType == 1
                error(strcat('Particle iP=', num2str(iP),' is too close to z=',num2str(b),' for DP one wall'));
            end
            iIm = iIm+1;
            imPts(iIm,:) = pts(iP,:); imPts(iIm,3) = 2*b - pts(iP,3);
            imForces(iIm,:) = 1.0*forces(iP,:);
            if has_torque && (b-pts(iP,3)) < dsep(2) && abs((b-pts(iP,3))-dsep(2)) > 1e-15
                iImt = iImt+1;
                imPtst(iImt,:) = imPts(iIm,:);
                imTorques(iImt,:) = 1.0*torques(iP,:);
            end
        end
    end
    
    % compute spread and interp matrices for monopoles
    [S,I,imS,imI] = SpreadInterpMatTrans(xE,yE,zpts,zwts,wm,phiMnorm,pts,imPts);
    % spread forces on particles
    [gridf,gridg,gridh] = SpreadForces(S,forces,imS,imForces,nxy,Nz);
    % get Fourier-Cheb coefficients of spread forces and derivatives
    [Cf,Cg,Ch,Df,Dg,Dh] = ftransformTransRHS(gridf,gridg,gridh,H);
    if has_torque
        % compute spread and interp matrices for dipole derivative
        [dSx,dIx,dSy,dIy,dSz,dIz,imdSx,imdIx,imdSy,imdIy,imdSz,imdIz] = ...
        SpreadInterpMatRot(xE,yE,zpts,zwts,wd,phiDnorm,phiDnorm_d,pts,imPtst);
        % spread curl of torques on particles
        [curlx,curly,curlz,dxT,dyT,dzT] = ...
          SpreadTorqueCurl(dSx,dSy,dSz,torques,imdSx,imdSy,imdSz,imTorques,nxy,Nz);
        % get Fourier-Cheb coefficients of spread curl and derivatives
        [CCurlf,CCurlg,CCurlh,DCurlf,DCurlg,DCurlh] = ftransformRotRHS(curlx,curly,curlz,H);
        % combine with trans RHS 
        Cf = Cf + CCurlf;
        Cg = Cg + CCurlg;
        Ch = Ch + CCurlh;
        Df = Df + DCurlf;
        Dg = Dg + DCurlg;
        Dh = Dh + DCurlh;
    end

    if nargin == 15
        if nargout < 14
            error('incorrect number of output variables requested');
        end
        varargout{1} = dIx; 
        varargout{2} = dIy; 
        varargout{3} = dIz; 
        varargout{4} = imdIx; 
        varargout{5} = imdIy; 
        varargout{6} = imdIz;
        varargout{7} = dxT;
        varargout{8} = dyT;
        varargout{9} = dzT;
        [St,It] = SpreadKernel(xE,yE,zpts,pts,zwts,wd,phiDnorm,2);
        [nIm,~] = size(imPtst);
        if nIm ~= 0
            [imSt,imIt] = SpreadKernel(xE,yE,zpts,imPtst,zwts,wd,phiDnorm,2);
        else
            imSt = []; imIt = [];
        end
        [Tx,~,~] = SpreadForces(St,torques,imSt,imTorques,nxy,Nz);
        varargout{10} = Tx;
        varargout{11} = It;
        varargout{12} = imIt;
    end
end
