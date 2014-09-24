function SHTransTest()

fileName = '../test/Original.out';

p = 6;
np = 2 * p * (p + 1);
fid = fopen(fileName,'r');
X = fscanf(fid,'%g');
fclose(fid);
X = reshape(X,3*np,[]);
X = transpose_lat_long(X,p);

subplot(2,2,1); plotb(X); title('original shape');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = '../test/ShtFiltering.out';

p = 6;
np = 2 * p * (p + 1);
fid = fopen(fileName,'r');
XX = fscanf(fid,'%g');
fclose(fid);
XX = reshape(XX,3*np,[]);
X  = transpose_lat_long(XX(:,1),p);
Y  = transpose_lat_long(XX(:,2),p);

subplot(2,2,2); plotb(X); title('pre filter');
subplot(2,2,3); plotb(Y); title('post filter');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = '../test/ShtResampling.out';

p = 2*6;
np = 2 * p * (p + 1);
fid = fopen(fileName,'r');
YY = fscanf(fid,'%g');
fclose(fid);
YY = transpose_lat_long(YY,p);
subplot(2,2,4); plotb(YY);title('resampled');
end

function X=transpose_lat_long(X,p)
X = reshape(X,[],3);
x = reshape(X(:,1),[],p+1)';
y = reshape(X(:,2),[],p+1)';
z = reshape(X(:,3),[],p+1)';
X = [x(:);y(:);z(:)];
end
