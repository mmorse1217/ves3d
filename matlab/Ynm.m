function [Yn Hn Wn] = Ynm(n, m, u, v)
% YNM(n,m,u,v) - returns the spherical harmonics of degree n and order m on
% the grid points u and v.
%
% n - degree
% m - order, if m is empty, it returns all orders (m=-n,...,n)
% u - north-south lines (0,pi)
% v - east-west lines   (0,2pi)
%
% SEE ALSO: parDomain.
%

  if(nargin==0), testYnm(); return; end
  np = length(u);
  if(isempty(m)), m =-n:n;end
  if(m>n), Yn = zeros(np,1); Hn = Yn; Wn = Yn; return; end

  Pn = legendre(n, cos(u))';
  nor = getNorConst(n);
  Pn = repmat(nor,length(u),1).*Pn(:,[n+1:-1:2 1:n+1]);

  %-- Spherical harmonics
  Yn = Pn.*exp(1i*v*(-n:n));
  Yn = Yn(:,n+1+m);
  if(nargout<2), return;end

  %-- First derivative
  Pn = [zeros(np,1) Pn zeros(np,1)];
  c = sqrt((n+(-n:n+1)).*(n+1-(-n:n+1)))/2;
  DPn = repmat(c(1:end-1),np,1).*Pn(:,1:end-2)...
	- repmat(c(2:end),np,1).*Pn(:,3:end);

  Hn = -DPn.*exp(1i*v*(-n:n));
  Hn = Hn(:,n+1+m);
  if(nargout<3), return;end

  %-- Second derivative
  DPn = [zeros(np,1) DPn zeros(np,1)];
  D2Pn = repmat(c(1:end-1),np,1).*DPn(:,1:end-2)...
         - repmat(c(2:end),np,1).*DPn(:,3:end);

  Wn = D2Pn.*exp(1i*v*(-n:n));
  Wn = Wn(:,n+1+m);


function nc = getNorConst(n)
% GETNORCONST - returns the normalizing constant for the Legendre
% polynomials of degree n
%
  persistent NOR

  if(n>=length(NOR) || isempty(NOR) || isempty(NOR{n+1}))
    NOR{n+1} = [(-1).^(-n:0) ones(1,n)].*sqrt((2*n+1)/4/pi *...
                                              myFact(n-abs(-n:n))./myFact(n+abs(-n:n)));
  end
  nc = NOR{n+1};


function val = myFact(N)
% MYFACT - tabulated factorial
%
  persistent factTable

  q = max(N);
  if(q>=length(factTable))
    tail = length(factTable);
    factTable = [factTable factorial(tail:q)];
  end
  val = factTable(N+1);

function testYnm()
  printMsg('Testing the spherical harmonic functions', 'sep', '-', 'indent', 1);
  p = 8;  
  [u v] = parDomain(p);
  [trash gwt]=grule(p+1);
  wt = pi/p*repmat(gwt', 2*p, 1);
  wt = wt(:);
  for l = 0:p
    Yl = Ynm(l,[],u,v);
    err = max(max(abs(Yl'*(Yl.*repmat(wt,1,2*l+1)) - eye(2*l+1))));
    printMsg('n=%d :%4.2e\n', l, err);
  end
  printMsg('','indent', -1);
  
