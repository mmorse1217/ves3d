clear all;clc

fileName = '../test/EvolveSurfMT_1.out';

p = 12;
nv = 1;
np = 2 * p * (p + 1);
fid = fopen(fileName,'r');
XX = fread(fid,'double');
fclose(fid);
XX = reshape(XX,3*np*nv,[]);

clf;
for ii=1:size(XX,2)
  PlotShape(reshape(XX(:,ii),[],nv),p);
  axis on;
  pause;
end

