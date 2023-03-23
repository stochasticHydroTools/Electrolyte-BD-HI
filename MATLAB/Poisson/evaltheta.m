% This function is used to evaluate the Chebyshev series at a given value
% of theta. 
% Inputs: phihat, the Chebyshev series (assumed to be a 3D array with
% Chebyshev coefficients along the third dimension), number of z points Nz,
% and angle to evaluate at, theta (somewhere between 0 and pi)
% Outputs: the value phi0 (a 2D matrix) of the Chebyshev series at that
% angle. 
function phi0 = evaltheta(phihat,Nz,theta)
    phi0=0*phihat(:,:,1);
    for j=0:Nz-1
        phi0=phi0+phihat(:,:,j+1)*cos(j*theta);
    end
end