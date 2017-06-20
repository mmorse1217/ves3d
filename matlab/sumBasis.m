function f = sumBasis(coeff, basisType, isReal, u, v)
  
  if(nargin==0), testThis(); return;end
  if(nargin<3), isReal = false; end
  
  p = sqrt(size(coeff,1))-1;
  nf = size(coeff, 2);
  
  if(nargin<4)
    [u v] = parDomain(p);
  end
  f = zeros(size(u,1),nf);
  np = size(f,1);
    
  switch basisType %Ynm, Hnm, Wnm
   case 'Ynm'
    for n=0:p
      Yn = Ynm(n,[],u,v);
      head = n^2;
      for m=1:2*n+1
        ind = head + m;
        f = f + repmat(coeff(ind,:),np,1).*repmat(Yn(:,m),1,nf);
      end
    end

   case 'Hnm'
    for n=0:p
      [Yn Hn]= Ynm(n,[],u,v);
      head = n^2;
      for m=1:2*n+1
        ind = head + m;
        f = f + repmat(coeff(ind,:),np,1).*repmat(Hn(:,m),1,nf);
      end
    end
    
   case 'Wnm'
    
    for n=0:p
      [Yn Hn Wn] = Ynm(n,[],u,v);
      head = n^2;
      for m=1:2*n+1
        ind = head + m;
        f = f + repmat(coeff(ind,:),np,1).*repmat(Wn(:,m),1,nf);
      end
    end

   otherwise 
    error('The basis type is not defined');
  end
  
  if(isReal), f = real(f); end
  
function testThis()
  
  for p=[6 12 24]
    shc = zeros((p+1)^2,1);
    shc(1:p^2) = rand(p^2,1);
    f = shSyn(shc);
    g = sumBasis(shc,'Ynm'); disp(max(abs(g-f)));
  end
