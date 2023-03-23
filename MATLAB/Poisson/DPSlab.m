% Doubly periodic Poisson inside of a dielectric slab
% The answers: phiCharges = average potential around each charge (phibar in
% report), ECharges = average electric field at each charge (Ebar in
% report)
addpath('../DPFunctions')
format long
% INPUT PARAMETERS
calcenergy=false;
epsratiot=1e-10;%'M';        % dielectric jump ratio at the top (epsilon_t=1/epsratiot). Use 'M' for metallic
epsratiob=1e-10;%'M';        % dielectric jump ratio at the bottom (epsilon_b=1/epsratiob)
gw=0.25;            % width of the Gaussian
H=75;              % height of the slab
L=100;                 % half the xy period (it goes from -L to L)
nxy = 32;           % choose grid size (Ewald param calculated from that)
nC=10;               % number of charges in the slab
rng(0);
pts=[2*L*(rand(nC,2)-0.5) (rand(nC,1)*(H-9*gw)+4.5*gw)]; % points where charges are
charges=zeros(nC,1); % vector of charge strengths
charges(1:2:end)= 1;
charges(2:2:end)= -1;
% rng(1);
% pts = [2*1*(rand(nC,2)-0.5)];
% pts(:,3)=[0.01;0.30; 1.0; 1.95];
% pts=[0 0 0.05+1e-3/2];
wallallimg=0; % whether to do all the images if doing a single wall (for tests)
gridfac=1.4; % grid spacing = g_t/gridfac
supp=12; % this is support in # of grid points (n_g in paper)
delta=1e-4;
if (nC > 10)
    gridfac=1.2;
    supp=10;
    delta=5e-4;
end
% Precomputations
if (epsratiob =='M' && epsratiot~='M')
    error('Exactly one metallic BC has to be on the TOP wall')
end
% Compute gtot from xy grid size
% p = 7;
% gtot = sqrt(1/(4*p^2)+gw^2);
% nxy = ceil(2*L/(gtot/gridfac));
% nxy = nxy+mod(nxy,2);
% hexy=2*L/nxy;
hexy = 2*L/nxy;
gtot = gridfac*hexy; % hexy = gtot/gridfac
p = 0.5*(gtot^2-gw^2)^(-1/2);
%Gaussian support
numSds = supp/2*hexy/gtot;       % how many standard deviations to truncate the Gaussian at
Hew = numSds*gtot; % 1/2 the distance where we include images
if (wallallimg) % include all images if testing single wall
    Hew = H/2;
end
% Near field truncation parameters
rposs=(0:1:100)*min([L H+4*gw])/100;
Gposs=GnearE(rposs,gw/sqrt(2),p)./GnearE(rposs,gw/sqrt(2),0); % new cutoff 
if (Gposs(end) > delta)
    error("Increase the Ewald parameter; interactions with more than 1 image")
else
    rnf= rposs(Gposs < delta);
    rnf=rnf(1); % pointwise near field cutoff
    rcut = rnf+numSds*gw; % pairwise near field cutoff
end
% z grid - spacing at middle is given by (6*Hew+H)/2*sin(pi/(Nz-1)). That
% has to match hxy
Nz = ceil(pi/2 *(H+6*Hew)/hexy);
sprintf("Doing Ewald splitting with %d xy pts and %d z pts, Hew=%f, rnf=%f",2*L/hexy,Nz, Hew,rnf)
% Get the free space answer if desired
% [phiFree, Efree] = slabFreeSpaceSol(pts,charges,gw,H,epsratiot,epsratiob);
% % Call the main solver routine to get the answer on the charges
DPSlabMain;
ECharges
% Check the answer by calling the solver without Ewald splitting
% Disabling this in the uploaded version because it doesn't work if gw is small
% DoublyPeriodicCharges; % for a check
% p1=phiCharges_NE-phi0_NE;
% E1=ECharges_NE;
% Xi0SolverWithCor;
% if (epsratiot=='M' || epsratiob=='M')
% else
% phiCharges_NE=phiCharges_NE-phi0_NE; % for dielectrics with no metal, phi(0)=0
% end
% erphi = phiCharges-phiCharges_NE;
% erE = (ECharges-ECharges_NE)/mean(sqrt(sum(ECharges_NE.*ECharges_NE,2)));
%[phiFree, Efree] = slabFreeSpaceSol(pts,charges,gw,H,epsratiot,epsratiob);
