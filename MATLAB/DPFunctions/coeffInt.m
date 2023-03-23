% Integrate the Chebyshev series to get the coefficients of the first
% integral. This function uses a loop. Equivalent formulation is to apply
% the firstIntegralMatrix. 
function fint = coeffInt(f,N)
    fint=zeros(N,1);
    fint(1)=f(N+2); % the c0 mode
    fint(2)=f(1)-0.5*f(3);
    for j=2:N-1
        fint(j+1)=1/(2*j)*(f(j)-f(j+2)*(j<N-1));
    end
end