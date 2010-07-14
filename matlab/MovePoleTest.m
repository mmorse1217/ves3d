addpath ../../Ves3DMat/src/
addpath ../../Ves3DMat/util/

%clear all;clc

% p = 12;
% if(~any(strcmp(who,'Rcell')) || length(Rcell) ~= p + 1)
%   [Rall Rcell] = GenerateShRotMats(p);
% end

% for lambda = linspace(0,2*pi,p)
% for m=1:p-1
%   D{m} = [cos(m * lambda) sin(m * lambda);...
%          -sin(m * lambda) cos(m * lambda)];
% end

% vRot = eye(2*p);
% vRot(2:2*p-1, 2:2*p-1) = blkdiag(D{:});
% vRot(2*p, 2*p) = cos(p*lambda);

% R = Rcell{4};
% X = boundary(p,'eightHump');
% X = reshape(X.cart.vecForm(),[],3);
% shc = shAnaReal(X);
% for ii=1:3
%   temp = reshape(shc(:,ii),p+1,[]);
%   temp = vRot * temp'; 
%   shcTemp = zeros(p+1,2*p);
%   for jj=0:p
%     lim = 2*jj+1;
%     if(jj==p), lim = 2*jj;end
%     shcTemp(jj+1,1:lim)= R{jj+1} * temp(1:lim,jj+1);
%   end
%   shcRot(:,ii) = shcTemp(:);
% end
% Xr = shSynReal(shcRot);
% plotb(Xr);
% pause(.1);
% end

fileName = '../test/MovePole.out';

p = 12;
np = 2 * p * (p + 1);
fid = fopen(fileName,'r');
XX = fread(fid, 3*np, 'float');
fclose(fid);

XX = reshape(XX,3*np,[]);

clf;
NV = 1;
for ii=1:size(XX,2)/NV
  PlotShape(XX(:,NV*(ii-1) + (1:NV)),p);
  pause(.1);
end

