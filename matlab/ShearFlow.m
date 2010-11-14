clear all;clc

fileName = '../scratch/ShearFlow_twoVes.txt';

p = 12;
nv = 2;
np = 2 * p * (p + 1);
XX = ReadAscii(fileName);

XX = reshape(XX,4*np*nv,[]);
tension = XX(3*np*nv+1:end,:);
XX = XX(1:3*np*nv,:);

for ii=1:size(XX,2)
  cla;
  disp(ii);
  PlotShape_Row(reshape(XX(:,ii),[],nv),p, reshape(tension(:,ii),[],nv));
  axis on;
  axis tight;
  pause(.1);
end

% fileName = '../test/EvolveSurf.out';

% p = 12;
% nv = 2;
% np = 2 * p * (p + 1);
% XX = ReadAscii(fileName);

% XX = reshape(XX,4*np*nv,[]);
% tension = XX(3*np*nv+1:end,:);
% XX = XX(1:3*np*nv,:);

% for ii=1:size(XX,2)
%   cla;
%   disp(ii);
%   PlotShape(reshape(XX(:,ii),[],nv),p, reshape(tension(:,ii),[],nv));
%   axis on;
%   pause(1);
% end