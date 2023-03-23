% This function computes the Fourier-Chebyshev coefficients of
% the spread forces on the grid as well as their derivatives in z

% Input:
%   gridf,gridg,gridh - spread forces on grid (output by SpreadForces.m)
%   H - half length of z grid (H = (b-a)/2)
% Output:
%   Cf,Cg,Ch - Fourier-Chebyshev coefficients of spread forces
%   Df,Dg,Dh - Fourier-Chebyshev coefficients of z derivative of spread forces 

function [Cf,Cg,Ch,Df,Dg,Dh] = ftransformTransRHS(gridf,gridg,gridh,H)
    [nxy,~,Nz] = size(gridf);
    U = zeros(nxy,nxy,2*Nz-2);
    Cf = ftransform(gridf,Nz,U);
    Cg = ftransform(gridg,Nz,U);
    Ch = ftransform(gridh,Nz,U);
    Df = chebCoeffDiff(Cf,Nz,H);
    Dg = chebCoeffDiff(Cg,Nz,H);
    Dh = chebCoeffDiff(Ch,Nz,H);
end
