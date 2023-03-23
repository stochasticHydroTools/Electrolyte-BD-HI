% ----------------------------------------------------------------------- %
% Author: Alex Kaiser & Bill Bao
% Date: Dec 2016
% Description: Flexible 5-pt Kernel for IBM
%
% Conditions:
% ---------------------------------------------------------------------
% | Support | Odd/Even | 1st Mom. | 2nd Mom. | 3rd Mom. | Sum Squares |
% ---------------------------------------------------------------------
% |  5-pt   |    -     |    x     |    x     |    x     |      x      |
% ---------------------------------------------------------------------
%
% Notes:
% (1) Letting K = 0 gives the so-called 'Old 6-point Kernel'.
% (2) Letting K = (38-sqrt(69))/60 gives the so-called 'New
% 5-point Kernel'. This choice of K ensures that the kernel has continuous
% first, second, and third derivatives. It is also the smallest value of K
% such that the kernel is everywhere non-negative.
% (3) This code is slow. However, the closed form of the kernel is really
% quite horrible. If speed becomes an issue, however, that will be the way
% to go, so it will be implemented eventually.
% ----------------------------------------------------------------------- %


function val = flex5pt(x,KK)

% The special value of K that gives a C^3 kernel
% KK = (38 - sqrt(69))/60; 

phi = @(r) (136 - 40*KK - 40*r.^2 + sqrt(2)*sqrt(3123 - 6840*KK + 3600*KK.^2 - 12440*r.^2 + 25680*KK*r.^2 - 12600*KK.^2*r.^2 + 8080*r.^4 - 8400*KK*r.^4 - 1400*r.^6))/280; 
  
% function is even, 
x = abs(x);
% compute r in [0,1]
r = (x < 0.5).*x + (x >= 0.5 & x < 1.5).*(x-1) + (x >= 1.5 & x < 2.5).*(x-2);
p = phi(r);
pm1 = (4 - 4*phi(r) - KK - 4*r + 3*KK*r - r.^2 + r.^3)/6;
pm2 = (-2 + 2*phi(r) + 2*KK + r - 3*KK*r + 2*r.^2 - r.^3)/12;
val = (x < 0.5).*p + (x >= 0.5 & x < 1.5).*pm1 + (x >= 1.5 & x < 2.5).*pm2;

% if abs(x) < 0.5
%     r = x; 
%     val = phi(r); 
% elseif abs(x) < 1.5
%     r = x - 1; 
%     val = (4 - 4*phi(r) - KK - 4*r + 3*KK*r - r.^2 + r.^3)/6;
% elseif abs(x) < 2.5
%     r = x - 2
%     val = (-2 + 2*phi(r) + 2*KK + r - 3*KK*r + 2*r.^2 - r.^3)/12; 
% else 
%     val = 0; 
% end 




