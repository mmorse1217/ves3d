function [f FmatRet] = filterSh(f,freq,useMat)
% FILTERSH(f,freq) - filters the function f such that the highest spherical
% harmonic frequency that is present after filtering is freq. Therefore,
% filterSh(f,p) where p is the current band limit of f does not modify f.
%

  persistent Fmat matFreqs

  mem = whos('Fmat');
  if(mem.bytes > 4e8), Fmat = []; matFreqs = []; end
  
  %%-- Checking the input
  if(nargin==0), testFilterSh(); return; end
  np = size(f,1); p = (sqrt(2*np+1)-1)/2;
  
  if(nargin<2 || isempty(freq)), freq = p;end
  if(fix(freq)~=freq), error('New frequency should be an integer'); end
  if(freq>p)
    error('The filter cutoff frequency should be smaller than current size');
  end
  if(nargin<3), useMat = true; end
  if(p>48), useMat = false; end

    %%-- Filtering
  FmatRet = [];
  if(useMat)
    if(size(matFreqs,2)>0)
      xIdx = find(matFreqs(1,:) == p);
      yIdx = find(matFreqs(2,xIdx) == freq);
    else
      xIdx = []; yIdx = [];
    end
    
    if(isempty(yIdx))
      Fmat{end+1} = getMat(p,freq);
      matFreqs(:,end+1) = [p;freq];
      xIdx = length(Fmat);
      yIdx = 1;
    end
    idx = xIdx(yIdx);
    ff = Fmat{idx}*f;
  else
    shc = shAna(f);
    shc((freq+1)^2+1:end,:) = 0;
    ff = shSyn(shc, false);
  end
  if( isreal(f)), ff = real(ff);end
  
function Fmat = getMat(p,freq)

  persistent queryUser

  np = 2*p*(p+1);
  printMsg('  * Reading filtering matrix for p=%d to q=%d:\n',p,freq,'indent',1);
  fileName = ['Fmat' num2str(p) '_'];
  [Fmat  rFlag] = readData(fileName,freq,[np np]);
  if(~rFlag)
    printMsg('  Failed\n');                              
    printMsg('  * Generating the filtering matrix for p=%d to q=%d.\n', ...
              p,freq,'indent',1);
     if(p>64)
      if(isempty(queryUser))
        queryUser = timeQuary([' Building the filter matrix for p=' num2str(p) ' may take ' ...
                            'a long time.\n Proceed with building? y/[n]?']);
        if(isempty(queryUser)), queryUser = 'n';end
      end
      if(strcmpi(queryUser,'n')), return; end
     end
    
    useMat = false; I = eye(np);
    Fmat  = filterSh(I, freq, useMat); writeData(fileName,freq, Fmat);
    printMsg('  Successful\n','indent',-2);
  else
    printMsg('  Successful\n','indent',-1);
  end
  
function testFilterSh()

  [u v] = parDomain(12);
  f  = Ynm(12,0,u,v) + Ynm(8,0,u,v);
  ff = filterSh(f,8);
  e = ff - Ynm(8,0,u,v);
  disp(max(max(abs(e))));

  f =  exp(sin(u).^4.*cos(u).*cos(4*v));
  ff = filterSh(f,12);
  e = ff - f;
  disp(max(max(abs(e))));

  f =  exp(sin(u).^4.*cos(u).*cos(4*v));
  ff1 = filterSh(f,7);
  ff2 = interpsh(interpsh(f,7),12);
  e = ff1 - ff2;
  disp(max(max(abs(e))));
