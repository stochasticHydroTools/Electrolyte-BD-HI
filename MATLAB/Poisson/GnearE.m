% Near field kernel for the electric field. 
function val=GnearE(r,gw,p)
    % Assuming epsilon=1
    val1=(-(exp(-(r.^2./(4*gw.^2)))./(gw*sqrt(pi))) + ...
        (2*exp(-(r.^2./(4*gw^2 + 1./p.^2))))./(sqrt(4*gw^2 + 1./p.^2)*sqrt(pi)))./...
        (4*pi.*r) - (-erf(r./(2*gw)) + erf(r./sqrt(4*gw^2 + 1./p.^2)))./(4*pi.*r.^2);
    % Taylor series for small r < gw^2
    val2 = (1./(24*gw.^3.*pi^(3/2)) - 1./(3*(4.*gw.^2 + 1./p.^2).^(3/2)*pi^(3/2))).*r + ...
        (-(1./(160*gw^5*pi^(3/2)))+ 1./(5*(4*gw^2 + 1./p.^2).^(5/2)*pi^(3/2))).*r.^3;
    val = val1.*(r>gw^2)+val2.*(r<=gw^2);
    val(isnan(val))=val2(isnan(val));
end

