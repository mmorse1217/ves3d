clear all;clc
addpath ~/school/projects/Ves3DMat/src/
addpath ~/school/projects/Ves3DMat/util/

%% Curvature flow
fileName = '../data/dumbbell_cart12';
p = 12;
fid = fopen(fileName,'r');
XX = fscanf(fid,'%g');
fclose(fid);

XX = reshape(XX,[],3);
X = d3Vec;
X.x = reshape(reshape(XX(:,1),2*p,[])',[],1);
X.y = reshape(reshape(XX(:,2),2*p,[])',[],1);
X.z = reshape(reshape(XX(:,3),2*p,[])',[],1);

S = vesicle(X);


% %% Curvature flow
% fileName = '../data/dumbbell_cart12';
% p = 12;
% fid = fopen(fileName,'r');
% XX = fscanf(fid,'%g');
% fclose(fid);
% 
% XX = reshape(XX,[],3);
% X = d3Vec;
% X.x = reshape(reshape(XX(:,1),2*p,[])',[],1);
% X.y = reshape(reshape(XX(:,2),2*p,[])',[],1);
% X.z = reshape(reshape(XX(:,3),2*p,[])',[],1);
% 
% S = vesicle(X);
% ts = .005;
% 
% for idx=1:30
%   S.cart = S.cart+ts*(S.geoProp.H.*S.geoProp.nor);
%   S.plot;axis equal;drawnow;
% end
%  
% fileName = 'cart12_final';
% p = 12;
% fid = fopen(fileName,'r');
% XX = fscanf(fid,'%g');
% fclose(fid);
% 
% XX = reshape(XX,[],3);
% X = d3Vec;
% X.x = reshape(reshape(XX(:,1),2*p,[])',[],1);
% X.y = reshape(reshape(XX(:,2),2*p,[])',[],1);
% X.z = reshape(reshape(XX(:,3),2*p,[])',[],1);
% 
% T = vesicle(X);
% figure;T.plot; axis equal;
% disp(max(norm(S.cart-T.cart)));