% solve poisson equation lap(phi) = -rho
% in triply periodic domain [-L,L]^3

% solver selection 
% (0 - spectral, 1 - 7pt Lap, 2 - 24pt Lap)
solver = 1;
% kernel selection
% (0 - gauss, 1 - ES 6pt, 2 - pesk 4pt, 3 - pesk 5pt, 4 - pesk 6pt)
kernType = 4;

% domain setup
N = 128;
L = 32;
h = 2*L/N;

% x and y grid
x = (-N/2:N/2-1)*h; y = x; z = x;
% wave nums
kvec = 2*pi*[0:N/2-1 -N/2:-1]/(2*L);
[Kx,Ky,Kz] = meshgrid(kvec);

% spectral
if solver == 0 
  D2 = Kx.^2 + Ky.^2 + Kz.^2;
  [Dx,Dy,Dz] = meshgrid(1i*kvec); 
% 7 pt lap Fourier symbol
elseif solver == 1 
  % Fourier symbol of second centered difference operator
  [D2x,D2y,D2z] = meshgrid((4/h^2)*sin(kvec*h/2).^2);
  D2 = D2x + D2y + D2z;
  % Fourier symbols of first centered difference operators
  [Dx,Dy,Dz] = meshgrid(1i.*sin(kvec.*h)/h);
% 24 pt lap Fourier symbol  
elseif solver == 2 
  lapl_24pt_FT = (4*cos(Kx*h)+4*cos(Ky*h)+4*cos(Kz*h) ...
                 + 4*cos(Ky*h).*cos(Kz*h)+4*cos(Kx*h).*cos(Kz*h) ...
                 + 4*cos(Kx*h).*cos(Ky*h)-24)/(6*h^2);
  D2 = lapl_24pt_FT./(1+h^2/12.*lapl_24pt_FT);
  [Dx,Dy,Dz] = meshgrid(1i.*(-(1/6)*sin(2*kvec*h)+4*sin(kvec*h)*(1/3))/h);
else
  error('solver enum invalid');
end 

% equal and opposite unit charges
charges = [1;-1];

if kernType == 0 % gaussian
    gw = 12*h/8;
    w = N;
    kern = @(x) 1/(sqrt(2*pi*gw^2))*exp(-x.^2/(2*gw^2));
elseif kernType == 1 % 6pt ES 
    w = 6;
    alpha = w*h/2;
    beta = 1.77*w;
    ker = @(z) ((z/alpha).^2 <= 1).*exp(beta*(sqrt(1-(z/alpha).^2)-1));%
    normM = integral(@(r)ker(r),-alpha,alpha);
    kern = @(z) ker(z)./normM;
elseif kernType == 2 % 4 pt Peskin
    w = 4;
    kern = @(r) 1/h*stnd4pt(r/h);
elseif kernType == 3 % 5pt Peskin
    w = 5; 
    kern = @(r) 1/h*flex5pt(r/h,(38-sqrt(69))/60);
elseif kernType == 4 % 6pt Peskin
    w = 6;
    kern = @(r) 1/h*flex6pt(r/h,59/60-sqrt(29)/20);
end

nR = 70;
nTrials = 25;
E = zeros(nR,nTrials);
nonC = E;

rng default;
for iR = 1:nR
    xrand1 = -L/8 + L/4 * rand;
    yrand1 = -L/8 + L/4 * rand;
    zrand1 = -L/8 + L/4 * rand;
    for iTrial = 1:nTrials
        utrans = randn(3,1); utrans = utrans./norm(utrans,2);
        xrand2 = xrand1 + utrans(1) * 0.1 * h * (iR-1); 
        yrand2 = yrand1 + utrans(2) * 0.1 * h * (iR-1); 
        zrand2 = zrand1 + utrans(3) * 0.1 * h * (iR-1); 
        pts = [xrand1,yrand1,zrand1;xrand2,yrand2,zrand2];
        R(iR) = norm(pts(1,:)-pts(2,:));
        Rhat = (pts(2,:)-pts(1,:))./R(iR);
        E1 = zeros(1,nTrials); nonC1 = E1;
        % spread charges onto grid
        [S,I] = SpreadWtsTP(x,y,z',pts,w,kern);
        rhogrid = permute(reshape(full(S*charges),N,N,N),[2,1,3]);
        % fourier transform rhs 
        rhohat = fftn(rhogrid);
        % solve poisson equation in fourier domain and ignore k=0
        phihat = -rhohat./D2; phihat(1,1,1) = 0;
        % compute electric field in Fourier space
        Exhat = -Dx.*phihat; 
        Eyhat = -Dy.*phihat;
        Ezhat = -Dz.*phihat;
        % back transform into real space
        phi = real(ifftn(phihat));
        Ex = real(ifftn(Exhat));
        Ey = real(ifftn(Eyhat));
        Ez = real(ifftn(Ezhat));
        % get field strength on charges
        Exq = interpolate(I,Ex,N,N);
        Eyq = interpolate(I,Ey,N,N);
        Ezq = interpolate(I,Ez,N,N);
        Evec = [Exq(1), Eyq(1), Ezq(1)];
        % get field magnitude
        E1(iTrial) = Evec*Rhat';
        nonC1(iTrial) = norm(Evec - E1(iTrial)*Rhat)/E1(iTrial);
        % get potential on charges
        potq = interpolate(I,phi,N,N);
    end
    E(iR,:) = E1;
    nonC(iR,:) = nonC1;
    disp(iR);
end
