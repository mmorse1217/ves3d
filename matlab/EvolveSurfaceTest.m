clear all;clc

% u = @(x,y) -cos(x).*sin(y);
% v = @(x,y) sin(x).*cos(y);

% [x y] = meshgrid(linspace(-5*pi,5*pi,100));
% uu = u(x,y);
% vv = v(x,y);

% streamslice(x,y,uu,vv,5);

fileName = '../examples/EvolveSurf.out';
%fileName = '../test/EvolveSurf.out';

p = 48;
nv = 1;
np = 2 * p * (p + 1);
fid = fopen(fileName,'r');
XX = fread(fid,'double');
fclose(fid);
XX = reshape(XX,3*np*nv,[]);

clf;
for ii=1:size(XX,2)
  disp(ii);
  PlotShape(reshape(XX(:,ii),[],nv),p);
  axis on;
  pause(1);
end