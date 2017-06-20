function dF = DmSH(F,dv,du,useMats)
% DMSH - Takes the derivative of function(s) defined on the sphere with
% respect to parameters through spherical harmonic transformation. Each
% column of F is a function defined on a Gauss-Legendre-uniform grid on the
% sphere. The grid may be produced by PARDOMAIN().
%
% DmSH(F,dv,du) : takes the derivative of the function F in the azimuth
% direction(east-west) dv times and in the elevation(south-north)
% direction du times (du<2).
% 
% DmSH(F,dv,du, useMats) : By default DmSH uses matrices for
% differentation, if they are not present, it generates and saves them for
% future use to the '/data/' folder. However, this may slow down the first
% evaluation of the function for a given size. Set useMats to false for
% faster evaluation for large vectors.
%
% SEE ALSO: PARDOMAIN, SHANA, SHSYN.
%

  persistent Du Duu Dv freq

  pMax = 49;
  mem = whos('Du');
  if(mem.bytes > 4e8), Du=[]; Duu=[]; Dv=[]; freqs = []; end

  %-- Checking the input
  if(nargin==0), testDmSH(); return;end
  if(du>2)
    error(['Only first and second order differentiation is permitted' ...
           ' for elevation']);
  end

  %-- Getting the size
  [d1 d2] = size(F);
  p = (sqrt(2*d1+1)-1)/2;
  if(p~=fix(p)),
    error(['The input function should be defined on a' ...
           ' Gauss-Legendre--uniform grid. Check the size of input.']);
  end
  if(nargin<4), useMats = true;end

  if(useMats && p<pMax)
    idx = find(freq==p);
    if(isempty(idx))
      [M1 M2 M3] = getMats(p);
      Du{end+1} = M1;
      Dv{end+1} = M2;
      Duu{end+1} = M3;
      freq(end+1) = p;
      idx = length(freq);
    end
    
    dF = F;
    for k=1:dv, dF = Dv{idx}*dF;end
    if(du==1)
      dF = Du{idx}*dF;
    elseif(du==2)
      dF = Duu{idx}*dF;
    end
  
  else
    printMsg('  * Direct differentiation via SHT for p=%d and (du,dv)=(%d,%d):\n',p,du,dv,'indent',1);
    [u v] = parDomain(p);
    SHC = shAna(F);
    if(dv>0)
      mMat = (1i*repmat(-p:p,p+1,1)).^dv;
      mMat = shrinkShVec(mMat(:));
      SHC = repmat(mMat,1,d2).*SHC;
    end

    switch du
     case 0
      dF = shSyn(SHC);
     case 1
      dF = sumBasis(SHC, 'Hnm');
     case 2 
      dF = sumBasis(SHC, 'Wnm');
    end
    printMsg('  Successful\n','indent',-1);
  end
  
  if( isreal( F ) ), dF = real( dF ); end
  
function [Du Dv Duu] = getMats(p)

  np = 2*p*(p+1);
  printMsg('  * Reading the differentiation matrices for p=%d:\n',p,'indent',1);
  [Du  m1] = readData('Du'  ,p,[np np]);
  [Dv  m2] = readData('Dv'  ,p,[np np]);
  [Duu m3] = readData('Duu' ,p,[np np]);
  
  if(m1 && m2 && m3)
    printMsg('  Successful\n','indent',-1);
    return;
  end

  printMsg('  Failed\n','indent',-1);
  printMsg('  * Generating the differentiation matrices for p=%d:\n', p,'indent',1);
  
  useMat = false; I = eye(np);
  if(~m1), Du  = DmSH(I, 0, 1, useMat); writeData('Du' ,p,Du );end
  if(~m2), Dv  = DmSH(I, 1, 0, useMat); writeData('Dv' ,p,Dv );end
  if(~m3), Duu = DmSH(I, 0, 2, useMat); writeData('Duu',p,Duu);end
  printMsg('  Successful\n','indent',-1);
  

function testDmSH()
% Self-test function for DmSH.
  for p=[8 12 24 32]
    [u v] = parDomain(p);
    f = exp(sin(u).^4.*cos(u).*cos(4*v));
    %-- reference
    r.fu  = sin(u).^3.*cos(4*v).*(5*cos(u).^2-1).*exp(sin(u).^4.*cos(u).* ...
                                                      cos(4*v));
    r.fv  = -4*sin(u).^4.*cos(u).*sin(4*v).*exp(sin(u).^4.*cos(u).*cos(4*v));
    r.fuu = sin(u).^2.*cos(4*v).*exp(sin(u).^4.*cos(u).*cos(4*v)).* ...
            (25*cos(u).^3-13*cos(u)+46*cos(u).^4.*cos(4.*v)- 60* ...
             cos(u).^6.*cos(4*v)+25*cos(u).^8.*cos(4*v)- 12*cos(u).^2.* ...
             cos(4*v)+cos(4*v));

    r.fuv = -4*sin(u).^3.*sin(4*v).*(5*cos(u).^2-1).*...
            exp(sin(u).^4.*cos(u).*cos(4*v)).*(1+sin(u).^4.*cos(u).*cos(4*v));
    r.fvv = 16*sin(u).^4.*cos(u).*exp(sin(u).^4.*cos(u).*cos(4*v)).* ...
            (-cos(4*v)+sin(u).^4.*cos(u).*sin(4*v).^2);

    %-- numrtical
    useMat = true;
    n.fu = DmSH(f,0,1,useMat);
    n.fv = DmSH(f,1,0,useMat);

    n.fuu = DmSH(f,0,2,useMat);
    n.fuv = DmSH(f,1,1,useMat);
    n.fvv = DmSH(f,2,0,useMat);

    %-- Comparison
    sep = repmat('-',1,60);
    flds = fieldnames(r);
    fprintf([sep '\n error for f = exp(sin(u)^4*cos(u)*cos(4v)), with n = %g\n' sep '\n'],p);
    for ii=1:length(flds)
      ff = flds{ii};
      err.(ff) = max(abs(r.(ff) - n.(ff)));
    end

    disp(err);
    disp(sep);
  end

  X = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
  plotb(X(:), f);
  title('plot of f over the sphere');
