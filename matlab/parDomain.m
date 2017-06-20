function [u v] = parDomain(p)
% PARDOMAIN - Returns the Gauss-Legendre--uniform grid on the unit sphere
% (0,pi)x(0,2pi).
%
% SEE ALSO: GRULE 
%
  
lambda = (0:2*p-1)'*pi/p;
theta = acos(grule(p+1));
[v u] = meshgrid(lambda,theta);
u = u(:); v = v(:);