% This is the main code for the doubly periodic solver with Gaussian cloud
% charges that DOES NOT USE EWALD SPLITTING
% Right now it is set up to compare directly to the DPSlab answer (using
% Ewald splitting). Uncomment as appropriate. 
% function [phiCharges_NE, ECharges_NE] = DoublyPeriodicCharges(nxy,supp,seed)
% addpath('../DPFunctions')
% % Define points (centers of Gaussian clouds) and charges
% gw=0.1; % width of the Gaussian
% % pts=[-0.1 -0.2 -0.3];
% % r0=[0.529316759100169,0.481570824029505,0.698507916905197];
% % % pts(2,:)=pts(1,:)+8*gw*r0;
% rng(seed);
% nC=100;
% pts = [2*(rand(nC,2)-0.5) rand(nC,1)*(2-9*gw)+4.5*gw]; % points where charges are
% charges = zeros(nC,1);
% charges(1:2:nC) = 1;
% charges(2:2:nC) = -1;
% epsratiot=1;
% epsratiob=1;
% Grid information
% The x and y grid
% L=1;            % grid goes from -L to L
nxy=128;        % number of xy points 
hexy=2*L/nxy;
xE=(-nxy/2:nxy/2-1)*hexy;
yE=(-nxy/2:nxy/2-1)*hexy;
[x,y]=meshgrid(xE,yE);
% Define the charge densities on the top and bottom as functions of x and y
sigmabot = -sum(charges)/2*exp(-(x.^2+y.^2)/(2*0.2^2))/(sqrt(4*pi^2*0.2^4));
sigmatop = -sum(charges)/2*exp(-(x.^2+y.^2)/(2*0.2^2))/(sqrt(4*pi^2*0.2^4));
% Fourier wave numbers in x and y
kvec=[0:nxy/2-1 -nxy/2:-1]/(2*L);
[kx,ky]=meshgrid(kvec,kvec);
ksq=(kx*2*pi).^2+(ky*2*pi).^2;
% The z grid
% H = 2;
Nz = ceil(pi/2*H/(hexy)); % match maximum spacing to hexy
[zpts,zwts] = clencurt(Nz-1);
zpts=H/2*(zpts+1);
zwts=H/2*zwts;
% Precompute the Chebyshev integrals on [-H,H]
[pints,cints]=precomputeInts(Nz,H/2);
% Distribute the charges on the grid by convolving with a Gaussian
supp=64;    % support of the Gaussian
[S,I] = SpreadWts(xE,yE,zpts,pts,supp,gw,zwts);
gridf = -spread(S,charges,nxy,Nz);
% First take the FFT3 to get frequency coefficients in x and y and
% Chebyshev coefficients in z
U=zeros(nxy,nxy,2*Nz-2);
fhats = ftransform(gridf,Nz,U);
sthat = fft2(sigmatop); % transforms of the wall charge densities
sbhat = fft2(sigmabot); % transforms of the wall charge densities
phihat=zeros(nxy,nxy,Nz);
SIMat = H^2/4*secondIntegralMatrix(Nz);
BCs = BCRows(Nz);
BCs(1,:) = BCs(1,:)*epsratiot;
BCs(3,:) = BCs(3,:)*epsratiob;
% Solve the BVPs for x and y with the proper BCs and RHS of the BCs
for iX=1:length(xE)
    for iY=1:length(yE)
        fhat=reshape(fhats(iY,iX,:),Nz,1);
        phihat(iY,iX,:)=SIMat(1:Nz,:)*BVPChebInt(sqrt(ksq(iY,iX)),...
            Nz,SIMat,BCs,H/2,fhat,sthat(iY,iX)*epsratiot*(ksq(iY,iX)~=0),...
            -sbhat(iY,iX)*epsratiob*(ksq(iY,iX)~=0));
    end
end
% For k=0 solveBVPSpectral solves the problem with homogeneous BCs
% Do the correction here by adding to the linear mode the integral fhat*z
% Complication because integrals cints and pints are on [-H/2,H/2], need to
% properly shift to get them on [0,H]
phihat(1,1,2) = phihat(1,1,2)- 0.5*(cints + H/2*pints)*reshape(fhats(1,1,:),Nz,1);
% Add the dipole moment of the charged walls
phihat(1,1,2) = phihat(1,1,2) + 0.5*H*sthat(1,1);
% Take the derivatives (E = -grad(phi)) in Fourier space
Exhat=-1i*2*pi*kx.*phihat;
Exhat(:,nxy/2+1,:)=0; % zero out the unpaired mode
Eyhat=-1i*2*pi*ky.*phihat;
Eyhat(nxy/2+1,:,:)=0; % zero out the unpaired mode
Ezhat = -chebCoeffDiff(phihat,Nz,H/2); % Chebyshev differentiation matrix for z
% Now IFFT3 back
phigrid = btransform(phihat,Nz,U);
Exgrid = btransform(Exhat,Nz,U);
Eygrid = btransform(Eyhat,Nz,U);
Ezgrid = btransform(Ezhat,Nz,U);
% Integrate back to get the potential on the charges
phi0_NE=phigrid(nxy/2+1,nxy/2+1,Nz); % evaluate pointwise potential at (0,0,0)
phiCharges_NE=interpolate(I,phigrid,nxy,Nz);
Exgrid=interpolate(I,Exgrid,nxy,Nz);
Eygrid=interpolate(I,Eygrid,nxy,Nz);
Ezgrid=interpolate(I,Ezgrid,nxy,Nz);
ECharges_NE=[Exgrid Eygrid Ezgrid];
% end