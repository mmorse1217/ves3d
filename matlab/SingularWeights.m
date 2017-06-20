function val  = SingularWeights(n)

  LP_of_Y = zeros(1, n+1);
  for k = 0:n
    LP_of_Y(k+1) = 4*pi./(2*k+1)*Ynm(k,0,0,0);
  end

  [trash gwt]=grule(n+1); 
  wt = pi/n*repmat(gwt',1,2*n); 
  wt = wt(:)';
  [u, v] = parDomain(n);
  Yf = 0; 
  for j = 0:n
    Yf = Yf + LP_of_Y(j+1)*(reshape(Ynm(j,0,u,v),2*n*(n+1),1)');
  end
  val = wt.*Yf./cos(u(:)'/2);