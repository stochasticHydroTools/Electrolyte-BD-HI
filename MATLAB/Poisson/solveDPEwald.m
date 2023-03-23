% Main solver for the doubly periodic problem within the spectral Ewald
% framework. 
% Inputs: f in real space on a grid, number of points in x and y (nxy),
% number of points in z (Nz), height of the slab H (this is a physical
% parameter), thet0=value of theta to get numbers at z =0 (so that phi mismatch
% at 0 is 0), and vector of wave numbers kvec. 
function [phihat, Exhat, Eyhat, Ezhat] = solveDPEwald(gridf,nxy,Nz,H,thet0,kvec,electroneutral)
    % First take the FFT3 to get frequency coefficients in x and y and
    % Chebyshev coefficients in z
    U=zeros(nxy,nxy,2*Nz-2);
    fhats = ftransform(gridf,Nz,U);
    phihat=zeros(nxy,nxy,Nz);
    SIMat = H^2*secondIntegralMatrix(Nz);
    BCs = BCRows(Nz);
    [kx,ky]=meshgrid(kvec,kvec);
    ksq=(kx*2*pi).^2+(ky*2*pi).^2;
    for iX=1:nxy
        for iY=1:nxy
            fhat=reshape(fhats(iY,iX,:),Nz,1);
            phihat(iY,iX,:)=SIMat(1:Nz,:)*...
                BVPChebInt(sqrt(ksq(iY,iX)),Nz,SIMat,BCs,H,fhat,0,0);
        end
    end
    % For k=0 solveBVPSpectral solves the problem with homogeneous BCs
    % We are putting off the correction until the correction solve - don't
    % do it here because system might not be charge neutral. If the system
    % is charge neutral, do the correction
    if (electroneutral==1)
        [pints,cints]=precomputeInts(Nz,H);
        A=(cints+H*pints)*reshape(fhats(1,1,:),Nz,1);
        phihat(1,1,2) = phihat(1,1,2)- 0.5*A;
    end
    % Subtract a constant so that the potential is 0 at 0 (NEW 06/03/2019)
%     phihat(1,1,1) = phihat(1,1,1) - evaltheta(phihat(1,1,:),Nz,thet0);
    % Take the derivatives in Fourier space (E= -grad(phi))
    Exhat=-1i*2*pi*kx.*phihat;
    Exhat(:,nxy/2+1,:)=0; % zero out the unpaired mode
    Eyhat=-1i*2*pi*ky.*phihat;
    Eyhat(nxy/2+1,:,:)=0; % zero out the unpaired mode
    Ezhat = -chebCoeffDiff(phihat,Nz,H); % Chebyshev differentiation matrix for z
end