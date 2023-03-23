% Function to get the harmonic correction on the inside of the slab.
% Inputs: vector of wave numbers, relevant z's. Epsilon ratios, four
% mismatch matrices, H, Hew (that is the width of the Gaussian cloud), and
% slope of the correction for the k=0 mode. 
% Outputs: correction (in Fourier/Chebyshev space) for phi, Ex, Ey, Ez. 
function [phihat, Exhat, Eyhat, Ezhat]= ...
correctslab(kvec,zpts,epsratiob,epsratiot,mph0,mEh0,mphH,mEhH,H,Hew,A0,B0) 
    nxy=length(kvec);
    Nz=length(zpts);
    [kx,ky]=meshgrid(kvec,kvec);
    ksq=(kx*2*pi).^2+(ky*2*pi).^2;
    phihat = zeros(nxy,nxy,Nz);
    Ezhat  = zeros(nxy,nxy,Nz);
    region = (zpts>-Hew & zpts<H+Hew); % the region we interpolate from
    for iX=1:nxy
        for iY=1:nxy
            k = sqrt(ksq(iY,iX));
            if (k > 0)
                phiIn=zeros(Nz,1);
                EIn=zeros(Nz,1);
                % Calculate the correction solution analytically inside the
                % region where we need it. 
                if (epsratiot=='M' && epsratiob=='M')
                    [phiIn(region),EIn(region)]=calcCorSolBothM(zpts(region),k,H,...
                    mph0(iY,iX),mphH(iY,iX));
                elseif (epsratiot=='M') % One wall 
                    [phiIn(region),EIn(region)]=calcCorSolMOne(zpts(region),k,H,...
                    mph0(iY,iX),mEh0(iY,iX),mphH(iY,iX),epsratiob);
                else
                    [phiIn(region),EIn(region)]=calcCorSol(zpts(region),k,H,...
                        mph0(iY,iX),mEh0(iY,iX),mphH(iY,iX),mEhH(iY,iX),...
                        epsratiot,epsratiob);
                end
                phihat(iY,iX,:)=ftransformz(phiIn,Nz);
                Ezhat(iY,iX,:)=ftransformz(-EIn,Nz);
            end
        end
    end
    phihat(1,1,:)=ftransformz((A0*zpts+B0).*region,Nz);
    Ezhat(1,1,:)=ftransformz(-A0*ones(Nz,1).*region,Nz);
    % Take the derivatives in Fourier space
    Exhat=-1i*2*pi*kx.*phihat;
    Exhat(:,nxy/2+1,:)=0; % zero out the unpaired mode
    Eyhat=-1i*2*pi*ky.*phihat;
    Eyhat(nxy/2+1,:,:)=0; % zero out the unpaired mode
end

function c = ftransformz(y,N)
    Vy = [y; flipud(y(2:N-1))];
    Uy = fft(Vy);
    c = [Uy(1); Uy(2:N-1)+Uy(end:-1:N+1); Uy(N)]/(2*N-2); % this is the cosine transform of y
    theta = pi*(0:(N-1))'/(N-1);
    er=max(abs(cos((0:N-1).*theta)*c-y));
    if (er > 1e-5)
        warning("Possible loss of accuracy in transforming correction coefficients")
    end
end

% Routine to calculate the correction from the mismatches. Based on the
% analytical solution and optimized for double precision arithmetic. 
function [inPhi, inE] = calcCorSol(z,k,H,mp0,mE0,mpH,mEH,epsratiot,epsratiob)
    eb = 1/epsratiob;
    et = 1/epsratiot;
%     inPhi=(-exp(k*z)*(-1 + et)*(mE0 - eb*k*mp0) + ...
%           exp((2*H-z)*k)*(1 + et)*(mE0 - eb*k*mp0) + exp((H-z)*k)*(-1 + eb)*...
%           (mEH + et*k*mpH) - exp(k*(H + z))*(1 + eb)*...
%           (mEH + et*k*mpH))/((-1 + eb + et - eb*et + exp(2*H*k)*...
%           (1 + eb)*(1 + et))*k);
    % This version corrected for underflow/overflow
    inPhi=(-exp(k*(z-2*H))*(-1 + et)*(mE0 - eb*k*mp0) + ...
          exp(-z*k)*(1 + et)*(mE0 - eb*k*mp0) + exp((-H-z)*k)*(-1 + eb)*...
          (mEH + et*k*mpH) - exp(k*(-H + z))*(1 + eb)*...
          (mEH + et*k*mpH))/((exp(-2*H*k)*(-1 + eb + et - eb*et) + ...
          (1 + eb)*(1 + et))*k);
%     inE = (-exp(k*z)*(-1 + et)*(mE0 - eb*k*mp0) - ...
%           exp((2*H-z)*k)*(1 + et)*(mE0 - eb*k*mp0) - exp((H-z)*k)*(-1 + eb)*...
%           (mEH + et*k*mpH) - exp(k*(H + z))*(1 + eb)*...
%           (mEH + et*k*mpH))/((-1 + eb + et - eb*et + exp(2*H*k)*...
%           (1 + eb)*(1 + et)));
    % This version corrected for underflow/overflow
    inE = (-exp(k*(z-2*H))*(-1 + et)*(mE0 - eb*k*mp0) - ...
          exp(-z*k)*(1 + et)*(mE0 - eb*k*mp0) - exp((-H-z)*k)*(-1 + eb)*...
          (mEH + et*k*mpH) - exp(k*(-H + z))*(1 + eb)*...
          (mEH + et*k*mpH))/(exp(-2*H*k)*(-1 + eb + et - eb*et) + ...
          (1 + eb)*(1 + et));
end

% Routine to calculate the correction from the mismatches for single
% metallic wall.
function [inPhi, inE] = calcCorSolMOne(z,k,H,mp0,mE0,mpH,epsratiob)
    epsb = 1/epsratiob; % and eps=1
    D = k*(exp(-k*H)*(1-epsb)+exp(k*H)*(1+epsb));
    Ai = (-(1+epsb)*k*mpH+exp(-k*H)*(-mE0+epsb*k*mp0))/D;
    Bi = ((-1+epsb)*k*mpH+exp(k*H)*(mE0-epsb*k*mp0))/D;
    inPhi = Ai*exp(k*z)+Bi*exp(-k*z);
    inE = k*Ai*exp(k*z)-k*Bi*exp(-k*z);
end

% Routine to calculate the correction from the mismatches for 2 metallic 
% boundaries. 
function [inPhi, inE] = calcCorSolBothM(z,k,H,mp0,mpH)
    Ai = (mpH-exp(-k*H)*mp0)/(exp(-k*H)-exp(k*H));
    Bi = (-mp0+mpH*exp(-k*H))/(1-exp(-2*k*H));
    inPhi = Ai*exp(k*z)+Bi*exp(-k*z);
    inE = k*Ai*exp(k*z)-k*Bi*exp(-k*z);
end