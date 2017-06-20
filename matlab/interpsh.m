function F = interpsh(F,newSize)
% Fint = interpsh(F,newSize), interpolates the set of functions F, defined
% on the sphere, to a grid of the new size newSize. The function is assumed
% to be defined on a Guess-Legendre--uniform grid, returned by PARDOMAIN.
% newSize is the number of spherical harmonics that is used to represent F.
% Therefore, newSize+1 will be the number of discretization points in the
% elevation direction. newSize can be larger or smaller than the current
% size of F. 
%
% SEE ALSO: PARDOMAIN, SHANA, SHSYN.
%

%-- self-test
  if(nargin==0), testInt(); return; end

  %-- Extracting sizes
  [d1 d2] = size(F);
  oldSize = (sqrt(2*d1+1)-1)/2;
  n = newSize; m = oldSize;

  if(newSize>oldSize)        %up sampling interpolation
    shc = zeros((n+1)^2, d2);
    shc(1:(m+1)^2,:) = shAna(F);
    F = shSyn(shc,isreal(F));
  elseif(newSize<oldSize)    %truncating the series
    shc = shAna(F);
    shc = shc(1:(n+1)^2,:);
    F = shSyn(shc,isreal(F));
  end
  
function testInt()

  n = 2; m = 1; p = 24;
  [u v] = parDomain(p);
  rho = real(Ynm(n,m,u,v));

  x = sin(u).*cos(v).*rho;
  y = sin(u).*sin(v).*rho;
  z = cos(u).*rho;

  subplot(1,3,2); plotb([x;y;z]); title('Original');

  X = interpsh([x y z],2*p); 
  subplot(1,3,3); plotb(X); title('Up sample');
  X = interpsh([x y z],p/2); 
  subplot(1,3,1); plotb(X); title('Down sample');
