function [qp qw] = grule(n)
% GRULE(n) - Returns the n-point Gauss-Legendre quadrature nodes and weights
%
% qp is the vector of quadrature points in [-1,1] interval and qw is the
% quadrature weight vector.
%
% This code is courtesy of L.N. Trefethen, Spectral method in Matlab, page 129.
%
% SEE ALSO: PARDOMAIN  
  
  if(nargin==0), testGrule();return;end
  
  beta = .5./sqrt(1-(2*(1:n-1)).^(-2));
  T = diag(beta,1) + diag(beta,-1);
  [V,D] = eig(T);
  qp = diag(D); 
  [qp,ind] = sort(qp); qp = qp';
  qw = 2*V(1,ind).^2;

function testGrule()
% TESTGRULE - Tester function for grule
%   
  m = 12;
  [qp qw] = grule(m);
  disp(['With ' num2str(m) ' points, polynomial of order...']);
  fprintf('\t -----------------------\n');
  for ii=[1:4:2*m-2 2*m-1:2*m+3]
    p = rand(1,ii+1);
    intN = sum(qw.*polyval(p,qp));
    
    p = [p./(ii+1:-1:1) 0];
    intA = polyval(p,1)-polyval(p,-1);
    fprintf('\t %-2g has error: %-2.3e\n',ii,abs(1-intN/intA));
	if(ii==2*m-1), fprintf('\t -----------------------\n');end
  end
  fprintf('\t -----------------------\n');