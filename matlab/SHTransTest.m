clear all;clc

fileName = '../test/ShtFiltering.out';

p = 12;
np = 2 * p * (p + 1);
fid = fopen(fileName,'r');
XX = fscanf(fid,'%g');
fclose(fid);
XX = reshape(XX,3*np,[]);

clf;
subplot(1,3,1); PlotShape(XX(:,1),p);
subplot(1,3,2); PlotShape(XX(:,2),p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = '../test/ShtResampling.out';

p = 24;
np = 2 * p * (p + 1);
fid = fopen(fileName,'r');
YY = fscanf(fid,'%g');
fclose(fid);
YY = reshape(YY,3*np,[]);

subplot(1,3,3); PlotShape(YY(:,1),p);

