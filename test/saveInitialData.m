% Script to make the Legendre and inverse Legendre transforms.

format long;
clear all;clc
addpath ../../Ves3DMat/src/ ../../Ves3DMat/util/

%% Options
precision = 'single';
p = 12;
shape = 'dumbbell';
filePath = '../data/';

actionList = {'InitialPos','legMat','rotMat'};

%% Save boundary
S = boundary(p,shape);

x = reshape(S.cart.x,p+1,[])';
y = reshape(S.cart.y,p+1,[])';
z = reshape(S.cart.z,p+1,[])';

X = [x(:); y(:); z(:)];
fileName = [filePath shape '_cart' num2str(p)];
saveData(fileName,X,precision,'bin');

%% Legendre matrix
[u w] = grule(p+1); u = u(:); u = acos(u);
v = zeros(size(u));

L  = cell(p+1,1);
LS = cell(p+1,1);
HS = cell(p+1,1);
WS = cell(p+1,1);

for n = 0:p
  [Yn Hn Wn] = Ynm(n, [], u, v);
  Yn = Yn*sqrt(2*pi);
  Hn = Hn*sqrt(2*pi);
  Wn = Wn*sqrt(2*pi);
	
  for m = 0:n
    L{1+m}(n+1,:) = w.*Yn(:,n+1+m).';
    LS{1+m}(:,n+1) = Yn(:,n+1+m);
    HS{1+m}(:,n+1) = Hn(:,n+1+m);
    WS{1+m}(:,n+1) = Wn(:,n+1+m);
  end
end

Lmat = []; LImat = []; Hmat = []; Wmat = [];
for m=0:p
  Lmat  = [Lmat ; reshape(L{m+1}(m+1:p+1,:),[],1)];
  LImat = [LImat; reshape(LS{m+1}(:,m+1:p+1),[],1)];
  Hmat  = [Hmat ; reshape(HS{m+1}(:,m+1:p+1),[],1)];
  Wmat  = [Wmat ; reshape(WS{m+1}(:,m+1:p+1),[],1)];
end

fileName = [filePath 'legTrans' num2str(p)]; 
saveData(fileName,Lmat,precision,'bin');
saveData(fileName,Lmat,precision,'txt');

fileName = [filePath 'legTransInv' num2str(p)]; 
saveData(fileName,LImat,precision,'bin');
saveData(fileName,LImat,precision,'txt');

fileName = [filePath 'd1legTrans' num2str(p)]; 
saveData(fileName,Hmat,precision,'bin');
saveData(fileName,Hmat,precision,'txt');

fileName = [filePath 'd2legTrans' num2str(p)]; 
saveData(fileName,Wmat,precision,'bin');
saveData(fileName,Wmat,precision,'txt');



% %% Test wrt the function shAna
% f = 1+6*cos(u); f = repmat(f,1,2*p);
% %-- DFT
% shc = fft(f,[],2)*sqrt(pi/2)/p;
% shc = shc(:,[p+1:end 1:p+1]); shc(:,[1 end]) = shc(:,[1 end])/2;
% for m =-p:p
% 	shc(:,p+1+m) = L{p+1+m}*shc(:,p+1+m);
% end
% shRef = reshape(shAna(f(:)),p+1,2*p+1);
% e = abs(shRef - shc);
% disp(max(max(e)));
% 
% %--inverse
% for m =-p:p
% 	ff(:,p+1+m) = LS{p+1+m}*shc(:,p+1+m);
% end
% ff(:,[1 end]) = 2*ff(:,[1 end]);
% ff = ff(:,[p+1:end 2:p]);
% ff = ifft(ff,[],2)*2*p/sqrt(2*pi);
% ff = real(ff);
% e = max(max(abs(f-ff)));
% disp(e);
