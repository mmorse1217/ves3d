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
    Ru = movePole(f,u(ii)-pi,0);
    Ru = Ru(ind,:);
    Rcell{1,ii} = Ru;
  end

  R = cell2mat(Rcell);