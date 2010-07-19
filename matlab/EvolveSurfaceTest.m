clear all;clc

fileName = '../test/EvolveSurf.out';

p = 12;
nv = 2;
np = 2 * p * (p + 1);
fid = fopen(fileName,'r');
XX = fread(fid,'float');
fclose(fid);
XX = reshape(XX,3*np*nv,[]);

clf;
for ii=1:size(XX,2)
  PlotShape(reshape(XX(:,ii),[],nv),p);
  pause(.1);
end

