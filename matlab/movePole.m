function shcRot = movePole(shc,thetaTarg,lambdaTarg) 
% [cc shcNew]= movePole(f,thetaTarg,lambdaTarg) move the pole of the
%parameters to the target point with coordinate thetaTarg (elevation
%[-pi/2,pi/2]) and lambdaTarg (azimuth [0,2pi)).  MOVEPOLE returns the
%Cartesian coordinates of the points with now parametrization (fRot) and
%their corresponding spherical harmonic coefficients (shcRot).

  
  if(nargin==0), testMovePole(); return;end
  N = size(shc,1)-1;
  
  alpha = lambdaTarg;
  beta = pi/2-thetaTarg; 
  % alpha rotation
  mMat = exp(1i*repmat(-N:N,N+1,1)*alpha);
  shc = mMat.*shc;

  % beta rotation
% Parametrization of the equator in the rotated from. We need 2*N+2
% equispaced points in [0,2pi).
  phiR = (0:2*N+1)*2*pi/(2*N+2)+pi/(2*N+2);
% Finding the corresponding material point
  [phi theta] = cart2sph(cos(beta)*cos(phiR),sin(phiR),-sin(beta)*cos(phiR));

  alphap = pi/3;   
%   phi = phi + pi;
% Mapping phi to [0,2pi) 
  phi = mod(phi,2*pi);

  % Coefficients for G
  c1 = -sin(beta)*sin(theta).*cos(phi) + cos(beta)*cos(theta);
  c2 = -sin(beta)*sin(phi)./cos(theta);

% Calculating f and g and finding their Fourier coefficients. 

  for n=0:N
    f = 0; g = 0;
    for m=-n:n
      Y = shBasis(n,m,phi,theta);       
      f = f + shc(n+1,N+1+m)*Y*exp(-1i*m*alphap);  % *** main difference ****%
      g = g + shc(n+1,N+1+m)*(c1.*dShBasis(n,m,phi,theta) + 1i*m*c2.*Y)*exp(-1i*m*alphap);
    end
%     ff(n+1,:) = fft(f);
%     gg(n+1,:) = fft(g);
    for m =-n:n
      ff(n+1,N+m+1) = sum(f.*exp(-1i*m*phiR))/(2*N+2);
      gg(n+1,N+m+1) = sum(g.*exp(-1i*m*phiR))/(2*N+2);
    end      
  end
  
  % Calculation new shc coefficients
  for n=0:N
    for m=-n:n
      P0 = shBasis(n,m,0,0);
      Q0 = dShBasis(n,m,0,0);
      shcRot(n+1,N+m+1) = (ff(n+1,N+m+1)*P0 + gg(n+1,N+m+1)*Q0)/(P0^2+Q0^2);
    end
  end
  
  % alpha rotation -- back
  mMat = exp(-1i*repmat(-N:N,N+1,1)*alpha);
  shcRot = mMat.*shcRot;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = shBasis(n, m, az, el)
  
  Pn = legendre(n, sin(el));
  nor = sqrt((2*n+1)/4/pi * factorial(n-abs(m))./factorial(n+abs(m)));
  if(m<0), nor = (-1)^m*nor;end
  
  if(n==0)
    Y = nor*Pn;
  else
    Y = nor*squeeze(Pn(abs(m)+1,:,:)).*exp(1i*m*az);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = dShBasis(n, m, az, el)  
  
  mo = m; m = abs(m);
  if(size(az,1)==1)
    Pn = zeros(n+2,size(az,2));
    Pn(1:n+1,:) = legendre(n,sin(el));
  else
    Pn = zeros(n+2,size(az,1),size(az,2));
    Pn(1:n+1,:,:) = legendre(n,sin(el));
  end
  
  if(m)  
    DPnm = (n+m)*(n-m+1)*squeeze(Pn(m,:,:))/2-squeeze(Pn(m+2,:,:))/2; 
  else
    DPnm = -squeeze(Pn(2,:,:)); 
  end
  nc = sqrt((2*n+1)/4/pi * factorial(n-m)/factorial(n+m));
  if(mo<0), nc = (-1)^m*nc;end
  H = nc*DPnm.*exp(1i*mo*az);
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function testMovePole()
  np = 64/4;
  
  [el az] = parDomainGauss(np);
  rho = @(theta,phi) 1+.1*cos(6*theta);%.4+.2*exp(10*real(Ynm(2,1,theta,phi)).^2);
  [x y z] = sph2cart(az,el,rho(el,az));

  u = pi/4; v = 0;
%   cx = shAna(x); cx = real(cx);
%   ccx = movePole(cx,u,v); ccx = real(ccx);
%   disp(cx./(ccx+(ccx==0)));
  
  xx = shSynA(movePole(shAnaA(x),u,v),1);
  yy = shSynA(movePole(shAnaA(y),u,v),1);
  zz = shSynA(movePole(shAnaA(z),u,v),1);
  
%   xx = shSyn(movePole(shAna(xx),u,v),1);
%   yy = shSyn(movePole(shAna(yy),u,v),1);
%   zz = shSyn(movePole(shAna(zz),u,v),1);

  x = x(:,[1:end 1]); xx = xx(:,[1:end 1]); 
  y = y(:,[1:end 1]); yy = yy(:,[1:end 1]);
  z = z(:,[1:end 1]); zz = zz(:,[1:end 1]);
  
%   x = x(1:np/2,1:np); xx = xx(1:np/2,1:np); 
%   y = y(1:np/2,1:np); yy = yy(1:np/2,1:np); 
%   z = z(1:np/2,1:np); zz = zz(1:np/2,1:np); 
  
  subplot(1,2,1); surf(x,y,z);   axis([-1 1 -1 1 -1 1]); 
  set(gca,'cameraposition',[10 10 8]); colormap summer;
  subplot(1,2,2); surf(xx,yy,zz);axis([-1 1 -1 1 -1 1]); 
  set(gca,'cameraposition',[10 10 8]); colormap summer;

  [phi theta r]=cart2sph(xx,yy,zz);
  e = abs(r - rho(theta,phi));
  disp(['error = ' num2str(max(e(:)))]);
  