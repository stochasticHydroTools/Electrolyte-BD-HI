% Near field kernel for phi. 
function val=Gnear(r,gw,p)
    % Assuming epsilon=1
    val1=(erf(r./(2.*gw))-erf(r./sqrt(4.*gw.^2+1./p.^2)))./(4*pi*r);
    % Taylor series for small r < gw^2
    val2 = ((1./(4*gw*pi^(3/2)))-1./(2*sqrt(4*gw^2+1./p.^2)*pi^(3/2))) + ...
        (-1./(48*gw^3*pi^(3/2))+1./(6*(4*gw^2+1./p.^2).^(3/2)*pi^(3/2)))*r.^2 + ...
        ((1./(640*gw^5*pi^(3/2))) - 1./(20*(4*gw^2 + 1./p.^2).^(5/2)*pi^(3/2)))*r.^4;
    val = val1.*(r>gw^2)+val2.*(r<=gw^2);
    val(isnan(val))=val2(isnan(val));
end

