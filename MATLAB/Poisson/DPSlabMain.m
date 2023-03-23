% Main code for the DP slab
% First solve is for z < 0 with the images
Hs = (H+6*Hew)/2; % This is half the domain
% Compute the grid on [-2Hew, H+2Hew]
nxy=2*L/hexy;
xE=(-nxy/2:nxy/2-1)*hexy;
yE=(-nxy/2:nxy/2-1)*hexy;
[x,y]=meshgrid(xE,yE);
sigwall=0.2;
sigmabot = -sum(charges)*exp(-(x.^2+y.^2)/(2*sigwall^2))/(sqrt(4*pi^2*sigwall^4));
sigmatop = -0*sum(charges)/2*exp(-(x.^2+y.^2)/(2*sigwall^2))/(sqrt(4*pi^2*sigwall^4));
kvec=[0:nxy/2-1 -nxy/2:-1]/(2*L);
[kx,ky]=meshgrid(kvec,kvec);
kvals=sqrt((2*pi*kx).^2+(2*pi*ky).^2);
% The z grid
[zpts,zwts] = clencurt(Nz-1);
zpts=Hs*(zpts+1)-3*Hew;
thet0=acos(3*Hew/Hs-1);
thetH=acos((3*Hew+H)/Hs-1);
zwts=Hs*zwts;

% Images for far field (within H_E of wall)
if (wallallimg==1)
    botpts = find(pts(:,3) < inf);
    toppts = find(pts(:,3) > -inf);
else
    toppts = find(pts(:,3) > H-2*Hew);
    botpts = find(pts(:,3) < 2*Hew);
end
timgpts = [pts(toppts,1:2) 2*H-pts(toppts,3)];
if (epsratiot=='M')
    timgchs=-charges(toppts);
else
    timgchs = -charges(toppts)*(1/epsratiot-1)/(1/epsratiot+1);
end
bimgpts = [pts(botpts,1:2) -pts(botpts,3)];
if (epsratiob=='M')
    bimgchs = -charges(botpts);
else
    bimgchs =-charges(botpts)*(1/epsratiob-1)/(1/epsratiob+1);
end  
C1 = unique([toppts; botpts]);
C2 = setdiff(1:nC,C1);

% Compute spreading matrix and spread weights
[S,I] = SpreadWts(xE,yE,zpts,[pts;timgpts; bimgpts],supp,sqrt(1/(4*p^2)+gw^2),zwts);
[Sew,Iew] = SpreadWts(xE,yE,zpts,pts,supp,sqrt(1/(4*p^2)),zwts);
gf_C1 = -spread(S(:,C1),charges(C1),nxy,Nz);
gf_C2 = -spread(S(:,C2),reshape(charges(C2),length(C2),1),nxy,Nz);
gf_C1star = -spread(S(:,nC+1:end),[timgchs; bimgchs],nxy,Nz);

% Guess solution for z < 0, z > H using C1 charges with image strength
[phihout, ~, ~, Ezhout] = solveDPEwald(gf_C1,nxy,Nz,Hs,thet0,kvec,wallallimg);

% Compute solution in slab using all charges and C1 images
[phihin, Exhin, Eyhin, Ezhin] = solveDPEwald(gf_C1+gf_C2+gf_C1star,...
    nxy,Nz,Hs,thet0,kvec,wallallimg);

if (epsratiot=='M' && epsratiob=='M') % metallic boundary
    phi_topSE = phi_BC_top*ones(nxy,nxy);% put the value here
    phi_tophat = fft2(phi_topSE);
    phi_botSE = phi_BC_bot*ones(nxy,nxy);
    phi_bothat = fft2(phi_botSE);
    mph0 = evaltheta(phihin,Nz,thet0)-phi_bothat;
    mEh0 = zeros(size(mph0));
    mphH = evaltheta(phihin,Nz,thetH)-phi_tophat;
    mEhH = zeros(size(mphH));
    % Linear mode
    B0 = -mph0(1,1);
    A0 = -(mphH(1,1)-mph0(1,1))/H;
elseif (epsratiot=='M') % one metallic boundary
    % Mismatch on top (phi only)
    phi_topSE = zeros(nxy,nxy);
    phi_tophat = fft2(phi_topSE);
    mphH = evaltheta(phihin,Nz,thetH)-phi_tophat;
    mEhH = zeros(size(mphH));
    % Mismatch on bottom (phi and E)
    mph0 = evaltheta(phihin,Nz,thet0)...
        -2/(1/epsratiob+1)*evaltheta(phihout,Nz,thet0);
    mEh0 = -evaltheta(Ezhin,Nz,thet0)...
        +2/(1/epsratiob+1)*evaltheta(Ezhout,Nz,thet0)/epsratiob+fft2(sigmabot);
    % Linear mode (initial solve has homogeneous Dirichlet BC)
    E=2/(1/epsratiob+1)*evaltheta(Ezhout(1,1,:),Nz,pi);
    A0=E/epsratiob-mEh0(1,1);
    B0 = -mphH(1,1)-A0*H;
else
    % Compute mismatch
    mph0 = evaltheta(phihin,Nz,thet0)...
        -2/(1/epsratiob+1)*evaltheta(phihout,Nz,thet0);
    mEh0 = -evaltheta(Ezhin,Nz,thet0)...
        +2/(1/epsratiob+1)*evaltheta(Ezhout,Nz,thet0)/epsratiob+fft2(sigmabot);
    mphH = evaltheta(phihin,Nz,thetH)...
        -2/(1/epsratiot+1)*evaltheta(phihout,Nz,thetH);
    mEhH = -evaltheta(Ezhin,Nz,thetH)...
        +2/(1/epsratiot+1)*evaltheta(Ezhout,Nz,thetH)/epsratiot-fft2(sigmatop);
    E=2/(1/epsratiob+1)*evaltheta(Ezhout(1,1,:),Nz,pi);
    C=2/(1/epsratiot+1)*evaltheta(Ezhout(1,1,:),Nz,0);
    A0=mean([E/epsratiob-mEh0(1,1) C/epsratiot-mEhH(1,1)]); % linear mode correction
    B0=0; % will be fixed later when potential is set to 0 at 0
end

% Do correction
mph0(kvals > pi/hexy) = 0;
mEh0(kvals > pi/hexy) = 0;
mphH(kvals > pi/hexy) = 0;
mEhH(kvals > pi/hexy) = 0;
[phich, Exch, Eych, Ezch]= correctslab(kvec,zpts,epsratiob,epsratiot,mph0,mEh0,mphH,mEhH,H,Hew,A0,B0);

% Back transform
U=zeros(nxy,nxy,2*Nz-2);
if (wallallimg==1)
    phigrid = btransform(phihin,Nz,U);
    Exgrid = btransform(Exhin,Nz,U);
    Eygrid = btransform(Eyhin,Nz,U);
    Ezgrid = btransform(Ezhin,Nz,U);
else
    phigrid = btransform(phich+phihin,Nz,U);
    Exgrid = btransform(Exch+Exhin,Nz,U);
    Eygrid = btransform(Eych+Eyhin,Nz,U);
    Ezgrid = btransform(Ezch+Ezhin,Nz,U);
end

% Interpolation
phiCharges = interpolate(I(:,1:nC),phigrid,nxy,Nz);
ECharges = [interpolate(I(:,1:nC),Exgrid,nxy,Nz) ...
    interpolate(I(:,1:nC),Eygrid,nxy,Nz)...
    interpolate(I(:,1:nC),Ezgrid,nxy,Nz)];
EFarPt = [interpolate(Iew(:,1:nC),Exgrid,nxy,Nz) ...
    interpolate(Iew(:,1:nC),Eygrid,nxy,Nz)...
    interpolate(Iew(:,1:nC),Ezgrid,nxy,Nz)];

% Add the near field, which includes C1, C1 images, and C2
% Images for near field (within rnf of slab)
topptsN = find(pts(:,3) > H-rnf);
botptsN = find(pts(:,3) < rnf);
timgptsN = [pts(topptsN,1:2) 2*H-pts(topptsN,3)];
if (epsratiot=='M')
    timgchsN=-charges(topptsN);
else
    timgchsN = -charges(topptsN)*(1/epsratiot-1)/(1/epsratiot+1);
end
bimgptsN = [pts(botptsN,1:2) -pts(botptsN,3)];
if (epsratiob=='M')
    bimgchsN = -charges(botptsN);
else
    bimgchsN = -charges(botptsN)*(1/epsratiob-1)/(1/epsratiob+1);
end  
nearpts = [pts; timgptsN; bimgptsN];
nearcharges = [charges; timgchsN; bimgchsN];
% Call near field function, truncating kernel at rcut = rnf+4*gw
rcut = 0;%RF###
[phiNear,ENear] = nearField(nearpts,nearcharges,nC,gw,p,2*L,rcut);
phiCharges=phiCharges+phiNear;
ECharges=ECharges+ENear;

% Pointwise potential at 0
% [~,I0] = SpreadWts(xE,yE,zpts,[0 0 0],supp,sqrt(1/(4*p^2)+(1e-10)^2),zwts);
[~,I0] = SpreadWts(xE,yE,zpts,[0 0 0; 0 0 Hs],supp,sqrt(1/(4*p^2)+(1e-10)^2),zwts);
phi0 = interpolate(I0,phigrid,nxy,Nz);
% phin0 = nearField([0 0 0; nearpts],[0;nearcharges],1,gw/sqrt(2),p,2*L,rcut);
phin0 = nearField([0 0 0; 0 0 H;nearpts],[0;0;nearcharges],1,gw/sqrt(2),p,2*L,rcut);
phi0 = phi0+phin0;
if (epsratiot=='M' || epsratiob=='M')
else
    % Subtract the potential at 0 if there are no metallic boundaries
%     phiCharges = phiCharges - phi0;
    phiCharges = phiCharges - phi0(1);
end

