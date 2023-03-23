% This file solves the 1D BVP using a Galerkin method. We compute the
% mobility matrix using the Galerkin solve. Also included are the functions
% and routines necessary to apply it fast. 
% This is all done over a Chebyshev grid. It is straightforward to
% use a uniform grid as an intermediary. 
% Set up the points and forces;
pts=[0.44 -0.35];
force=[1;-1];
hg = 0.10; % width of the Gaussian

% Parameters 
N=32;   % number of Chebyshev points
k=2;
H=2;

% % Precomputations - the basis functions and nodes
n=0:N-1;
% a=(k*H + n.^2)./(4 + k*H + 4*n + n.^2); % values of a for the 0 RHS
a=zeros(1,N);
[y,w]=clencurt(N-1);
[yup,wup]=clencurt(2*N+2); % has to be one less for clencurt call
E=cos(acos(y).*n)-a.*cos(acos(y).*(n+2));
Eup=cos(acos(yup).*n)-a.*cos(acos(yup).*(n+2));
Ehat = spdiags([ones(N,1) -a'], [0 -2], N,N);
Ehat(N-1,N-1)=1-a(N-1); Ehat(N-2,N)=-a(N);
Ehat_up =  spdiags([ones(N,1) -a'; zeros(N+3,2)], [0 -2], 2*N+3,2*N+3);

% % Precompute the matrices A and B
DEUp = zeros(2*N+3,N);
for iCol=1:N
    coeffs = reshape(full(Ehat_up(:,iCol)),1,1,2*N+3);
    d1cos = chebCoeffDiff(coeffs,2*N+3,H);
    DEUp(:,iCol) = btransformz(reshape(d1cos,2*N+3,1),2*N+3);
end
Anum = DEUp'*diag(wup*H)*DEUp+...
    k*(Eup(1,:)'*Eup(1,:)+Eup(2*N+3,:)'*Eup(2*N+3,:));
Bnum = Eup'*diag(wup*H)*Eup;
D = Anum+k^2*Bnum;

% All other matrices for the order N^2 method
R = evalup(N,2*N+3);    % upsampling matrix
Psi=R'*diag(H*wup)*R;
S = exp(-(y-pts).*(y-pts)/(2*hg^2))/sqrt(2*pi*hg^2);    % spread matrix to Chebyshev nodes
J = S'; % interpolation matrix (adjoint of S)
Inv_Gal = E*D^(-1)*E'*Psi;
MGalerkin = J*Psi*Inv_Gal*S;    % FULL MOBILITY MATRIX
pVel = MGalerkin*-force;

% % This part is about applying M fast (in linear time)
Psif = applyRT(diag(wup*H)*applyR(S*-force,N),N);
Etpsif = applyET(ftransformz(Psif,N),N,Ehat);
sol = applyE(D^(-1)*Etpsif,N,Ehat);
pVel_FAST = J*applyRT(diag(wup*H)*applyR(sol,N),N);

% Leslie's solver - compute mobility column by column
MLeslie = zeros(length(pts));
for iC=1:2
    force=zeros(2,1);
    force(iC)=1;
    fhat = ftransformz(S*-force,N);
    SI = H^2*secondIntegralMatrix(N);
    BCs = BCRows(N);
    uhat = SI*BVPChebInt(k,N,SI,BCs,H,fhat,0,0);
    u = btransformz(uhat,N);
    MLeslie(iC,:) = J*diag(w*H)*u;
end

% Test problem with inhomogeneous BC - this is just a BVP with right hand
% side "RHS." This does NOT involve spread and interpolate. 
% RHS = (k^2+1)*sin(y*H);
% gplus = cos(H)+k*sin(H);
% gminus = cos(-H)-k*sin(-H);
% Pf = applyRT(diag(wup*H)*applyR(RHS,N),N);
% EtPf = applyET(ftransformz(Pf,N),N,Ehat);
% EtPf = EtPf+applyET(ftransformz([gplus; zeros(N-2,1); -gminus],N),N,Ehat);
% ugal = applyE(D^(-1)*EtPf,N,Ehat);
% utr = sin(y*H);
% er=abs((utr-ugal));


% FUNCTIONS TO COMPUTE THE MATRICES / APPLY THEM FAST

% Compute the dense upsampling matrix R
function R = evalup(Nin,Nout)
    theta_in = pi*(0:(Nin-1))'/(Nin-1);
    inmat = (cos((0:Nin-1).*theta_in));
    theta_out = pi*(0:(Nout-1))'/(Nout-1);
    outmat = (cos((0:Nin-1).*theta_out));
    R = outmat*(inmat^(-1));
end

% Fast way of applying R
function fup = applyR(f,N)
    fhat = ftransformz(f,N);
    % Pad with zeros
    fhat(N+1:2*N+3) = 0;
    % Backtransform
    fup = btransformz(fhat,2*N+3);
end

% Fast way of applying R'
function fdown = applyRT(fup,N)
    fhat = btransformz(fup,2*N+3);
    % Truncate
    fhat(N+1:2*N+3)=[];
    % Back transform
    fdown = ftransformz(fhat,N);
end

% Fast way of applying E
function out = applyE(in,N,Ehat)
    outhat = Ehat*in;
    out = btransformz(outhat,N);
end

% Fast way of applying E'
function pIPs = applyET(fhat,N,Ehat)
    oddsums = sum(fhat(1:2:end));
    evensums = sum(fhat(2:2:end));
    IP=zeros(N,1);
    IP(1:2:N) = oddsums+(N-1)/2*fhat(1:2:end);
    IP(2:2:N) = evensums+(N-1)/2*fhat(2:2:end);
    IP(1) = IP(1) + (N-1)/2*fhat(1);
    IP(N) = IP(N) + (N-1)/2*fhat(N);
    % Now assign the components of IP 
    pIPs = Ehat'*IP;
end

% Go to Chebyshev space from real space
function c = ftransformz(y,N)
    Vy = [y; flipud(y(2:N-1))];
    Uy = fft(Vy);
    c = [Uy(1); Uy(2:N-1)+Uy(end:-1:N+1); Uy(N)]/(2*N-2); 
    % this is the cosine transform of y
end

% Real space back to Chebyshev space
function fgrid = btransformz(fhat,Nz)
    U=zeros(2*Nz-2,1);
    U(1)=fhat(1);
    U(2:Nz-1)=fhat(2:Nz-1)/2;
    U(Nz)=fhat(Nz);
    U(Nz+1:2*Nz-2)=fhat(Nz-1:-1:2)/2;
    U=ifft(U)*(2*Nz-2);
    fgrid = U(1:Nz);
end
