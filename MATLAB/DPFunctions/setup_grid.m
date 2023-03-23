% Function to set up the grid given x-y box length L and particle radius Rh
% for several kernel types. This works for no wall, bottom wall and 
% slit channel geometries.

% Input: L  - length in x,y
%        Rh - particle radius
%        h_min - min particle height (or approximate loc of bottom bndry) 
%        h_max - max particle height (or approximate loc of top bndry)
%        domType - 0 for no wall, 1 for bottom wall, 2 for slit channel
%                - if no wall, z \in [h_min-1.5*wh/2,h_max+1.5*wh/2] (see msep in routine)
%                - if one wall, z \in [0,h_max+1.5*wh/2]
%                - if two wall, z \in [h_min,h_max]
%        kernType - 0 for Gaussian (truncated at (approx) 5*sigma from mean)
%                 - 1 for ES 6 pt (both monopole and dipole if has_torque = true) 
%                 - 2 for ES 5 pt (both monopole and dipole if has_torque = true)
%                 - 3 for ES 4 pt (only monopole, dipole not supported) 
%                 - 4 for Peskin 6 pt
%                 - 5 for Peskin 5 pt
%                 - 6 for Peskin 4 pt
%                 - 7 for ES 6 pt monopole, ES 5 pt dipole
%                 - 8 for ES 6 pt monopole, ES 4 pt dipole
%                 - 9 for ES 5 pt monopole, ES 4 pt dipole 
%        varargin{1} = has_torque - boolean to specify whether problem involves torque
%                                 - NOTE: can't use 3 < kernType < 7 with has_torque=true


% Output:  xE   - periodic x grid [-L,L-hxy]
%          yE   - periodic y grid [-L,L-hxy]
%          zpts - Chebyshev grid [0,2H] (in reverse)
%          hxy  - grid spacing in x,y to match Rh
%          nxy  - length of xE, yE
%          Nz   - length of zpts
%          H    - half length of z domain
%          L    - corrected length of x,y domain (for integer grid size)
%          w    - number of points underneath kernel
%               - if has_torque = false, w = wm
%               - if has_torque = true, w = [wm,wd]
%          msep = [msep_down,msep_up] - max separation b/w monopole and lower/upper bndry
%               - if domType = 0, msep_down = msep_up = 1.5*wh/2
%               - if domType = 1, msep_down = wh/2, msep_up = 1.5*wh/2 
%               - if domType = 2, msep_down = msep_up = wh/2
%               - when images are applicable, if a particle is within msep of lower/upper
%                 wall, a monopole image will be used (in SpreadAndTransformSimple.m)
%          varargout{1} = dsep - max separation b/w dipole and lower/upper bndry

function [xE,yE,zpts,zwts,hxy,nxy,Nz,H,L,w,msep,varargout] = setup_grid(L,Rh,h_min,h_max,domType,kernType,varargin)
  if nargin > 6
    has_torque = varargin{1};
    if kernType > 3 && kernType < 7 && has_torque
      error('cannot create both monopole and dipole Peskin kernels (yet)')
    end
    if has_torque && nargout ~= 12
      error('incorrect number of outputs with torque');
    end
  else
    has_torque = false;
    if nargout ~= 11
      error('incorrect number of outputs without torque');
    end
  end
  % translation-only
  if ~has_torque
    switch kernType
      case 0 % gaussian w/o truncation
        gw = Rh/sqrt(pi); % stdev
        hxy = gw*2/3; % 12 grid pts per 8 stdev
        w = 15; % (truncate at 5*sigma = w*hxy/2 = w*sigma/3 => w = 15)
      case 1 % ES 6 pt
        w = 6;
        hxy = Rh/1.5539; % read from table 2
      case 2 % ES 5 pt
        w = 5;
        hxy = Rh/1.3437; % read from table 2
      case 3 % ES 4 pt
        w = 4;
        hxy = Rh/1.2047; % read from table 2
      case 4 % Peskin 6 pt
        w = 6;
        hxy = Rh/1.56; % read from table 4
      case 5 % Peskin 5 pt 
        w = 5;
        hxy = Rh/1.26; % read from table 4
      case 6 % Peskin 4 pt
        w = 4;
        hxy = Rh/1.1453; % read from table 4
      otherwise
        error('translation-only kernel type not supported');
    end
  % translation and rotation
  else
    switch kernType
      case 0 % gaussian
        gw = Rh/sqrt(pi); % stdev of monopole
        gwt = Rh/(6*sqrt(pi))^(1./3.); % stdev of dipole
        hxy = gwt*2/3; % 12 grid pts per 8 stdev of skinner Gaussian
        wd = 15; % (10*sigma_d = w*hxy = w*sigma_d*2/3 => w = 15)
        % width for monopole corresponds truncation at ~5.1*gw from mean 
        wm = 19; % (10*sigma_f = w*hxy = w*sigma_d*2/3 = w*sigma_f*sigma_d/sigma_f => w ~ 19)
        w = [wm,wd];
      case 1 % ES 6 pt monopole and dipole
        w = [6,6];
        hxy = Rh/1.7305; % read from table 3
      case 2 % ES 5 pt monopole and dipole
        w = [5,5];
        hxy = Rh/1.5598; % read from table 3
      case 7 % ES 6 pt monopole, 5 pt dipole
        w = [6,5];
        hxy = Rh/1.7309; % read from table 3
      case 8 % ES 6 pt monopole, 4 pt dipole
        w = [6,4];
        hxy = Rh/1.6121; % read from table 3
      case 9 % ES 5 pt monopole, 4 pt dipole
        w = [5,4];
        hxy = Rh/1.5382; % read from table 3
      otherwise
        error('translation-rotation kernel type not supported');
    end
  end

  % number of points in x,y
  nxy = 2*L/hxy;
  % correct L to give integer number of points
  if floor(nxy) ~= nxy
    fprintf('Adjusting requested L for integer grid size..\n');
    nxy = floor(nxy);
    L = hxy/2*nxy;
  end
  
  % half support with safety factor
  msep_down = 1.5*w(1)*hxy/2*(domType==0) + w(1)*hxy/2*(domType~=0);
  msep_up = 1.5*w(1)*hxy/2*(domType~=2) + w(1)*hxy/2*(domType==2);
  msep = [msep_down,msep_up];

  if has_torque
    dsep_down = 1.5*w(2)*hxy/2*(domType==0) + w(2)*hxy/2*(domType~=0);
    dsep_up = 1.5*w(2)*hxy/2*(domType~=2) + w(2)*hxy/2*(domType==2);
    dsep = [dsep_down,dsep_up];
    varargout{1} = dsep;
  end

  % left/right endpts for z grid and scaling factors for [a,b] -> [-1,1]
  a = h_min - msep_down*(domType==0); b = h_max + msep_up*(domType~=2);
  bma = b-a; bpa = b+a; H = bma/2;
  % upsample factor for cheb grid
  NzFac = 1.2;
  % number of points in z to match mid spacing with hxy with some upsampling
  Nz = ceil(NzFac*pi/(acos(-hxy/H)-pi/2));
  % x and y grid
  xE = (-nxy/2:nxy/2-1)*hxy;
  yE = (-nxy/2:nxy/2-1)*hxy;
  % chebyshev nodes and weights
  [zpts,zwts] = clencurt(Nz-1);
  % rescaling nodes and weights
  zpts = H*zpts + bpa/2;
  zwts = H*zwts;

  fprintf('Final grid settings for ');
  if domType == 0
    fprintf('No Wall domain:\n');
  elseif domType == 1
    fprintf('Bottom Wall domain:\n');
  elseif domType == 2
    fprintf('Two Wall domain:\n');
  end
  fprintf('\t N_xy = %d\n',nxy);
  fprintf('\t h_xy = %d\n',hxy);
  fprintf('\t L = %.16f\n',L);
  fprintf('\t N_z = %d\n',Nz);
  fprintf('\t z_a = %.16f \n \t z_b = %.16f\n',a,b);
  
  % determine L to give nearest power of 2 or 3
  nxyPow2_u = 2^(ceil(log2(nxy)));
  nxyPow2_d = 2^(floor(log2(nxy)));
  nxyPow3_u = 3^(ceil(log(nxy)/log(3)));
  nxyPow3_d = 3^(floor(log(nxy)/log(3)));
  Lpow2_u = nxyPow2_u*hxy/2;
  Lpow2_d = nxyPow2_d*hxy/2;
  Lpow3_u = nxyPow3_u*hxy/2;
  Lpow3_d = nxyPow3_d*hxy/2;
  if (nxy ~= nxyPow2_u && nxy ~= nxyPow2_d && ...
      nxy ~= nxyPow3_u && nxy ~= nxyPow3_d)
    fprintf('Try rerunning with L = %f so that N_xy = %d is a power of 2\n',Lpow2_d,nxyPow2_d);
    fprintf('Try rerunning with L = %f so that N_xy = %d is a power of 2\n',Lpow2_u,nxyPow2_u);
    fprintf('Try rerunning with L = %f so that N_xy = %d is a power of 3\n',Lpow3_d,nxyPow3_d);
    fprintf('Try rerunning with L = %f so that N_xy = %d is a power of 3\n',Lpow3_u,nxyPow3_u);
  end
end
