function [SHC DFT LT]= shAnaReal(F)
  
%-- Extracting the size
  [d1 d2] = size(F);
  p = (sqrt(2*d1+1)-1)/2;
  if(p~=fix(p))
    error(['The input should be defined on a Gauss-Legendre--uniform'...
           ' grid. Check the size of input.']);
  end
  
  DFT = getDFT(p);
  LT = getLegMat(p);
  
  for sur=1:d2
    shc = zeros(p+1,2*p);
    f = reshape(F(:,sur),p+1,2*p)';   
    f = DFT * f;
    f = f';
    
    for m = 0:p
      ind = 2*m:2*m+1;
      if(m == 0), ind = 1;end
      if(m == p), ind = 2*p;end
      
      shc(:,ind) = LT{m+1} * f(:,ind);
    end
    SHC(:,sur) = shc(:);
  end
  
function LT = getLegMat(p)

  [el wt]=grule(p+1); 
  u = acos(el(:));

  LT  = cell(p+1,1);
  for n = 0:p
    LT{n+1} = zeros(p + 1);
  end
  
  for n = 0:p
    Yn = Ynm(n, [], u, 0*u) * sqrt(2*pi); %at v=0, Ynm = Pnm;
    for m = 0:n
      LT{1+m}(n+1,:) = wt.*Yn(:,n+1+m).';
    end
  end
  
function MAT = getDFT(p)

  MAT = zeros(2*p);
  grid = (0:2*p-1)*pi/p;
   
  MAT(:,1) = 1/2;
  for jj = 1:p   
    MAT(:,2*jj) = cos(jj*grid);
    if( jj == p), break; end
    MAT(:,2*jj + 1) = sin(jj*grid);    
  end
  MAT(:,end) = MAT(:,end)/2;
  MAT = MAT' * grid(2)/pi;

  