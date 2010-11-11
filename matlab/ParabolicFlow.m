clear all;clc

fileName = '../scratch/ParabolicFlowII.out';
%fileName = '../scratch/ParabolicFlowOffset_u01_ts05.txt';

p = 12;
nv = 1;
np = 2 * p * (p + 1);
XX = ReadAscii(fileName);

XX = reshape(XX,4*np*nv,[]);
tension = XX(3*np*nv+1:end,:);
XX = XX(1:3*np*nv,:);

for ii=1:size(XX,2)
  cla;
  x = reshape(XX(:,ii),[],3);
  x = x(:,1);
  rr(ii,:) = [min(x) max(x) max(x)-min(x)];
  disp([ii min(x) max(x) max(x)-min(x)]);
  PlotShape(reshape(XX(:,ii),[],nv),p, reshape(tension(:,ii),[],nv));
  axis on;
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