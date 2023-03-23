% This function computes the Fourier-Chebyshev coefficients of
% the spread curl of torques on the grid as well as their derivatives in z

% Input:
%   curlx,curly,curlz - spread curl of torques on grid (output by SpreadTorqueCurl.m)
%   H - half length of z grid (H = (b-a)/2)
% Output:
%   CCurlf,CCurlg,CCurlh - Fourier-Chebyshev coefficients of spread curl of torques
%   DCurlf,DCurlg,DCurlh - Fourier-Chebyshev coefficients of z derivative of spread curl of torques 

function [CCurlf,CCurlg,CCurlh,DCurlf,DCurlg,DCurlh] = ftransformRotRHS(curlx,curly,curlz,H)
    [nxy,~,Nz] = size(curlx);
    U = zeros(nxy,nxy,2*Nz-2);
    CCurlf = 1/2*ftransform(curlx,Nz,U);
    CCurlg = 1/2*ftransform(curly,Nz,U);
    CCurlh = 1/2*ftransform(curlz,Nz,U);
    DCurlf = chebCoeffDiff(CCurlf,Nz,H);
    DCurlg = chebCoeffDiff(CCurlg,Nz,H);
    DCurlh = chebCoeffDiff(CCurlh,Nz,H);
end
