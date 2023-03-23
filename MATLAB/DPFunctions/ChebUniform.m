% ChebUniform.m - this is a script which showcases going back and forth
% between the Chebyshev and uniform grids on [-H, H]
%% Initialization
Nc = 128;    % number of Chebyshev points
Nu = 128;    % number of uniform points
H=2;        % height
[xc,wc]=clencurt(Nc-1); % Chebyshev nodes and weights
xc=xc*H;    % Transform to [-H, H]
wc=wc*H;
hu=2*H/(Nu-1); % uniform point spacing
xu = (-H:hu:H)';
% The function we want to sample (a Gaussan)
syms x
h=0.1;
f = 1/sqrt(2*pi*(2*h)^2)*exp(-x^2/(2*(2*h)^2));
% Values on the Chebyshev grid
fc = double(subs(f,xc));
% Values on the uniform grid
funif =double(subs(f,xu));
W=diag(wc);

%% Chebyshev to uniform
% Compute the matrix I
% First compute the ftransform matrix G
C2U = chebtoUni(Nc,xc/H,xu/H);
% Use C2U to evaluate on the uniform grid
fu = C2U*fc;
% Compute the error
erI1=abs(fu-funif);
% % Now the other way - trying to use the adjoint but this never worked
% S1 = W^(-1)*I1'*h;
% % Check the error in spreading
% erS1=abs(S1*funif - fc);

%% Uniform to Chebyshev
% Extend the uniform grid
supp = 16; % number of points for interpolant
xuu=(-H-4*hu:hu:H+4*hu)';
funif =double(subs(f,xuu));
U2C = unitoCheb(xc,xuu,supp,hu);
% Use U2C to evaluate on the Chebyshev grid and compute the error. 
erS2=abs(U2C*funif-fc);
% % Now compute the interpolation as the adjoint of spreading - this never
% worked. 
% I2 = S2'*W/h;
% erI2=abs(I2*fc-funif);


%% Functions to form the dense matrices

% From Chebyshev to uniform. Note that the inputs must be sets of points on
% [-1, 1]. 
function I1 = chebtoUni(Nc,xc,xu)
    G = cos((0:Nc-1).*acos(xc))^(-1);
    B = cos((0:Nc-1).*acos(xu));
    I1 = B*G;
end

% Matrix from the extended uniform grid to the Chebyshev grid. Note that
% the input xu is the EXTENDED uniform grid. 
function S2 = unitoCheb(xc,xu,supp,h)
    inds=1:length(xu);
    Nc = length(xc);
    S2=zeros(length(xc),length(xu));
    for iC=1:Nc % loop over the Chebyshev points (rows of the spreading matrix)
        % Find the points near
        iInds = inds(abs(xu-xc(iC)) < supp/2*h);
        % Build the Lagrange interpolant
        for iInd=1:length(iInds)
            j = iInds(iInd);
            S2(iC,j) = prod(xc(iC)-xu(setdiff(iInds,j)))./prod(xu(j)-xu(setdiff(iInds,j)));
        end
    end
end