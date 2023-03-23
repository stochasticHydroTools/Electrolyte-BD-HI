% Precompute the integrals for the Nz Chebyshev modes that come up in the
% A_p A_x and A_y values 
function [pints, uvints] = precomputeInts(Nz,H)
    theta = pi*(0:1000-1)'/(1000-1);
    [~,zwts2]= clencurt(999);
    zwts2=H*zwts2;
    pints=zwts2*(cos((0:Nz-1).*theta)); % T_j(z) 
    pints=pints(1:Nz);
    uvints=H*zwts2*(cos((0:Nz-1).*theta).*cos(theta)); % z* T_j(z)
    uvints=uvints(1:Nz);
end