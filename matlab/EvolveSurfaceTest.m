clear all;clc

fileName = '../test/EvolveSurf.out';

p = 12;
np = 2 * p * (p + 1);
fid = fopen(fileName,'r');
XX = fscanf(fid,'%g');
fclose(fid);
XX = reshape(XX,3*np,[]);

clf;
NV = 1;
for ii=1:size(XX,2)/NV
  PlotShape(XX(:,NV*(ii-1) + (1:NV)),p);
  pause(.1);
end

