function f = shSyn(shc,isReal)
% F = SHSYN(SHC,isReal) - Calculates the inverse spherical transform for the
% given set of spherical harmonics coefficients SHC. Each column of SHC
% is the coefficients of a function. isReal is the flag that indicates
% that the desired function F is real valued.
%
% SEE ALSO: SHANA, YNM.
%

  persistent LTI LTI_Freqs

  mem = whos('LTI');
  if(mem.bytes > 4e8), LTI = []; LTI_Freqs = []; end
  
  %-- Checking the input
  if(nargin==0), testShSyn(); return; end
  if(nargin==1), isReal = false; end

  shc = expandShVec(shc);
  %-- Calculating the sizes
  [d1 d2]= size(shc);
  p = (sqrt(8*d1+1)-3)/4;

  if(p~=fix(p))
    error(['The input should be defined on a Gauss-Legendre--uniform'...
           ' grid. Check the size of input.']);
  end

  %-- Inverse Legendre transform
  idx = find(LTI_Freqs == p);

  if(isempty(idx))
    LTI{end+1} = getLegSynMat(p);
    LTI_Freqs(end+1) = p;
    idx = length(LTI);
  end

  f = zeros(2*p*(p+1),d2);
  for m = -p:p-1 %The last frequency is extra
    ind = (m+p)*(p+1)+1;
    f(ind:ind+p,:) = LTI{idx}{p+1+m}*shc(ind:ind+p,:);
  end

  f(1:p+1,:) = 2*f(1:p+1,:);
  f = circshift(f,-p*(p+1));
  f = reshape(f,p+1,2*p,d2);
  f = ifft(f,[],2)*2*p;
  if(isReal), f = real(f);end
  f = reshape(f,[],d2);

function LTI = getLegSynMat(p)

  printMsg('  * Reading the inverse Legendre Transform matrix for p=%d:\n',p,'indent',1);
  [LTI  rFlag] = readData('invLegTrans',p,[p+1 (2*p+1)*(p+1)]);
  
  if(~rFlag)
    printMsg('  Failed\n');
    printMsg('  * Generating the inverse Legendre Transform matrix for p=%d:',p,'indent',1);
    
    [el wt]=grule(p+1); 
    u = acos(el(:));

    LTI  = cell(1,2*p+1);
    for n = 0:p
      Yn = Ynm(n, [], u, 0*u); %at v=0, Ynm = Pnm;
      for m = -n:n
        LTI{p+1+m}(:,n+1) = Yn(:,n+1+m);
      end
    end
    writeData('invLegTrans',p, cell2mat(LTI));
    printMsg('  Successful\n','indent',-1);
  else
    LTI = mat2cell(LTI, p+1, repmat(p+1,1,2*p+1));
    printMsg('  Successful\n','indent',-1);
  end

function testShSyn()

  N = 8;
  [el az] = parDomain(N);

  sep = '-------------------------------\n';
  fprintf(['        ' sep '\t   n        m     error\n        ' sep]);
  for n=0:N
    for m=-n:n
      f = Ynm(n,m,el,az);
      s = shAna(f);
      ff = shSyn(shAna(f));
      fprintf('\t   %d        %d     %2.4e\n',n,m,max(max(abs(f-ff))));
    end
  end
  fprintf(['        ' sep ' \t A more general function \n        ' sep]);
  f = 1+cos(6*el);
  ff = shSyn(shAna(f),1);
  fprintf(['\t   error = %2.4e\n        ' sep],max(max(abs(f-ff))));
