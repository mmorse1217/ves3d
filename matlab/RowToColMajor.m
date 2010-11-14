function xr = RowToColMajor(xv,p)

  for ii=1:size(xv,2)
    X = reshape(xv(:,ii),[],3);
    x = reshape(X(:,1),2*p, p+1)';
    y = reshape(X(:,2),2*p, p+1)';
    z = reshape(X(:,3),2*p, p+1)';
    
    xr(:,ii) = [x(:);y(:);z(:)];
  end