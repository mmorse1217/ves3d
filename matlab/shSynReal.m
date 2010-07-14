function [F IDFT ILT]= shSynReal(SHC)
  
  if(nargin == 0) testThis();return; end;
  
  %-- Extracting the size
  [d1 d2] = size(SHC);
  p = (sqrt(2*d1+1)-1)/2;
  if(p~=fix(p))
    error(['The input should be defined on a Gauss-Legendre--uniform'...
           ' grid. Check the size of input.']);
  end
  
  %-- Transform
  ILT = getInvLegMat(p);
  IDFT = getIDFT(p);
    
  for sur=1:d2
    f = zeros(p+1,2*p);
    shc = reshape(SHC(:,sur),p+1,2*p);
       
    for m = 0:p
      ind = 2*m:2*m+1;
      if(m == 0), ind = 1;end
      if(m == p), ind = 2*p;end
    
      f(:,ind) = ILT{m+1} * shc(:,ind);
    end
    
    %-- DFT
    f = f';
    f = IDFT * f;
    f = f';
    F(:,sur) = f(:);
  end

function LT = getInvLegMat(p)

  [el wt]=grule(p+1); 
  u = acos(el(:));

  LT  = cell(p+1,1);
  for n = 0:p
    LT{n+1} = zeros(p + 1);
  end
  
  for n = 0:p
    Yn = Ynm(n, [], u, 0*u)* sqrt(2*pi); %at v=0, Ynm = Pnm;
    for m = 0:n
      LT{1+m}(n+1,:) = Yn(:,n+1+m).';
    end
  end
  
  for m = 0:n
      LT{1+m} = LT{1+m}';
  end
  
function MAT = getIDFT(p)

  MAT = zeros(2*p);
  grid = (0:2*p-1)*pi/p;
   
  MAT(:,1) = 1;
  for jj = 1:p   
    MAT(:,2*jj) = cos(jj*grid);
    if( jj == p), break; end
    MAT(:,2*jj + 1) = sin(jj*grid);    
  end

function testThis()  
  p = 12;
  
  for ii=1:p+1
    for jj=1:2*p
      shc = zeros(p+1,2*p);     
            
      shc(ii,jj) = 1;             
      f = shSynReal(shc(:));
      c2 = shAnaReal(f);
      disp(max(abs(c2-shc(:))));
    end
  end