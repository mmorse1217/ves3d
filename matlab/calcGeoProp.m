function [gp gpAll]= calcGeoProp(S)
% gp = calcGeoProp(X)  calculates the geometric properties of the srurface
% defined by X. X is of class vesicle. gp is a structure containing the
% first fundamental coefficients E,F,G, second fundamental coefficients
% L,M,N, surface normal vector nor, area element W, mean and Gaussian
% curvatures H,K, spherical Laplacian Lap, as well as the function handle
% for surface Gradient and surface divergence.
%

  % NOTES: Bellow, gpAll holds the geoProps in the upsampled frequency. The
  % fundamental forms are not continuous and cannot be interpolated via
  % spherical harmonics, it is more accurate to calculate them in the original
  % shOrder. The geoProps in the same spherical order of S (or X) are stored in
  % gp.
  %
  % Also note that interpolating the normal causes it to lose its unity size
  % (n.n=1); if you decide to interpolate, you need to renormalize.
  %
  if(nargin==0), testGP(); return; end

  %-- Function handles for filtering and upsampling depending on the options
  % fH is the filtering handle and 
  % dH is the differentiation handle 
  % X is upsampled based on options
  if(S.upFreq>S.p && S.filterFreq==S.p)
    X = interpsh(S.cart, S.upFreq); 
    dH = @(f,dv,du) DmSH(interpsh(f, S.upFreq), dv, du);
    fH = @(f) interpsh(f, S.p);
  elseif(S.upFreq>S.p)
    X = interpsh(S.cart, S.upFreq); 
    dH = @(f,dv,du) DmSH(interpsh(f, S.upFreq), dv, du);
    fH = @(f) filterSh(interpsh(f,S.p), S.filterFreq);
  else
    X = S.cart;
    dH = @DmSH;
    fH = @(f) filterSh(f, S.filterFreq);
  end
  
  %-- Taking derivatives (in pUp)
  gpAll.Xu  = X.DmSH(0,1); gpAll.Xv  = X.DmSH(1,0);
  gpAll.Xuu = X.DmSH(0,2); gpAll.Xuv = X.DmSH(1,1); gpAll.Xvv = X.DmSH(2,0);

  %-- First fundamental coeffs
  gpAll.E = dot(gpAll.Xu, gpAll.Xu); 
  gpAll.F = dot(gpAll.Xu, gpAll.Xv); 
  gpAll.G = dot(gpAll.Xv, gpAll.Xv);
  W2 = abs(gpAll.E.*gpAll.G - gpAll.F.^2); gpAll.W = sqrt(W2);
  gpAll.nor = cross(gpAll.Xu,gpAll.Xv); gpAll.nor = gpAll.nor./gpAll.W;
  
  %-- Second fundamental coefficients, not generally continuous so cannot be
  %interpolated down
  gpAll.L = dot(gpAll.Xuu, gpAll.nor); 
  gpAll.M = dot(gpAll.Xuv, gpAll.nor);
  gpAll.N = dot(gpAll.Xvv, gpAll.nor);
  
  %-- Gauss and mean curvatures
  gpAll.H = 0.5*(gpAll.E.*gpAll.N - 2*gpAll.F.*gpAll.M + gpAll.G.*gpAll.L)./W2;
  gpAll.K = (gpAll.L.*gpAll.N - gpAll.M.*gpAll.M)./W2;
    
  %-- Mapping back to p 
  gp.H = fH(gpAll.H); gp.K = fH(gpAll.K);

  %-- Calculating the discontinuous functions in p frequency
  gp.Xu = S.cart.DmSH(0,1);
  gp.Xv = S.cart.DmSH(1,0);
  % normal is continuous, but interpolating it will cause it to loose its unit
  % length. Either we can calculate in p or renormalize after downsampling. We
  % calculate it in p, because it is Xu and Xv are spectrally accurate.
  gp.nor = cross(gp.Xu, gp.Xv); 
  gp.W = norm(gp.nor);
  gp.nor = gp.nor./gp.W;
  
  %-- Surface operators coefficients
  Cu = (gpAll.G.*gpAll.Xu - gpAll.F.*gpAll.Xv)./W2;
  Cv = (gpAll.E.*gpAll.Xv - gpAll.F.*gpAll.Xu)./W2;

  gp.Grad_noisy = @(f) dH(f,0,1).*Cu + dH(f,1,0).*Cv;
  gp.Div_noisy  = @(F) sum(dH(F,0,1).*Cu + dH(F,1,0).*Cv,2);

  %-- Surface operators
  gp.Grad = @(f) fH(gp.Grad_noisy(f));
  gp.Div  = @(F) fH(gp.Div_noisy(F));
  gp.Laplacian  = @(f) fH(gp.Div_noisy(gp.Grad_noisy(f)));

  %-- Linearized curvature
  gp.linCurv = @(X) fH(( gpAll.E.*dot(dH(X,2,0),gpAll.nor) - ...
                       2*gpAll.F.*dot(dH(X,1,1),gpAll.nor) + ...
                         gpAll.G.*dot(dH(X,0,2),gpAll.nor))./W2/2);

  %-- Linearized bending and tension operators on the surface
  H2K = 2*(gpAll.H.^2-gpAll.K);
  gp.bendingCoeff = @(f) -S.kappa*(gp.Laplacian(f) + fH(interpsh(f, S.upFreq).*H2K));
  gp.bendingOp = @(f) gp.bendingCoeff(f).*gp.nor;
  gp.tensionOp = @(f) gp.Grad(f) + 2*fH(interpsh(f,S.upFreq).*gpAll.H).*gp.nor;
%   gp.tensionOp = @(f) fH(gp.Grad_noisy(f) + interpsh(f,S.upFreq).*LapX);

  %Ordering fields 
  gp = orderfields(gp);
  
function testGP()

  ll = d3Vec;
  sep = [repmat('-',1,50) '\n'];
  printString = sprintf([sep '    p   \\int(H n)\t\\int(K)-4pi\t  Div(n) \n' sep]);
    
  for p = [4 8 12 16 24]

    X = boundary(p,'neck');
    X.upFreqCoeff = 2;
    X.cart = X.movePole(p/2,p/2);
    %- int(Hn)
    e1 = norm(integrateOverS(X, X.geoProp.H .* X.geoProp.nor));
    %- int(K)
    e2 = abs(integrateOverS(X, X.geoProp.K) - 4*pi);
    %- Div(n)
    HH = -X.geoProp.Div(X.geoProp.nor)/2;
    e3 = max(abs(HH-X.geoProp.H))/max(abs(X.geoProp.H));
    printString = [printString sprintf('   %d \t %4.2e \t %4.2e\t %4.2e\n',p,e1,e2,e3)];
  end
  fprintf('%s', printString);
