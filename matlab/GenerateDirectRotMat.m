function [R Rcell] = GenerateDirectRotMat(p)

  addpath ../../Ves3DMat/src/
  addpath ../../Ves3DMat/util/

  np = 2*p*(p+1);
  ind = reshape(1:np,p+1,[])';
  ind = ind(:);
  f = eye(np);
  f = f(:,ind);
  u = parDomain(p);
  Rcell = [];
  
  for ii=1:p+1
    disp(['Generating idx = ' num2str(ii)]);
    Ru = movePole(f,u(ii)-pi,0);
    Ru = Ru(ind,:);
    Rcell{1,ii} = Ru;
  end

  if ( nargout == 0 )
    fprintf('\n Testing :\n\n');
    X = boundary(p,'dumbbell');
    X = reshape(X.cart.vecForm(),[],3);
    X = X(ind,:);

    for ii=1:p+1
      disp(['target idx = ' num2str(ii)]);
      Xr(ind,:) = Rcell{1,ii} * X;
      plotb(Xr);
      pause;
    end
  else
    R = cell2mat(Rcell);
  end
  