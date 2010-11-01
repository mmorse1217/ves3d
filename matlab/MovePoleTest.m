function [] = MovePoleTest(p, task,precision)

  addpath ../../Ves3DMat/src/
  addpath ../../Ves3DMat/util/

  if ( nargin < 2 )
    task = 'DataCheck';
  end
  
  switch task
   case 'DataCheck'
    fileName = '../test/MovePole.out';

    nv = 1;
    np = 2 * p * (p + 1);
    fid = fopen(fileName,'r');

    XX = fread(fid,precision);
    fclose(fid);

    XX = reshape(XX,3*np*nv,[]);
    for jj=1:size(XX,2)
      clf;
      disp(jj);
      PlotShape(reshape(XX(:,jj),[],nv),p);
      drawnow;
      pause(.1);
    end

   case 'CheckOperators'
    [trash Rsh] = GenerateShRotMats(p);
    [trash Rdi] = GenerateDirectRotMat(p);
    
    np = 2 * p * ( p + 1);
    ind = reshape(1:np,p+1,[])';
    X = boundary(p,'dumbbell');
    X = reshape(X.cart.vecForm(),[],3);
    shc = shAnaReal(X);
    X_rowmajor = X(ind,:);

    for latitude = 1:p+1
      %- sh rotation
      R = Rsh{latitude};

      for ii=1:size(shc,2)
        temp = reshape(shc(:,ii),p+1,[]).';
        shcTemp = zeros(p+1,2*p);
        for jj=0:p
          lim = 2*jj+1;
          if(jj==p), lim = 2*jj;end
          shcTemp(jj+1,1:lim)= R{jj+1} * temp(1:lim,jj+1);
        end
        shcRot(:,ii) = shcTemp(:);
      end
      Xr_shc = shSynReal(shcRot);
      
      Xr_shc = Xr_shc(ind,:);
      subplot(1,2,1); PlotShape(Xr_shc(:),p);

      %- Direct rotation
      R = Rdi{latitude};
      Xr_direct = R*X_rowmajor;
    
      subplot(1,2,2); PlotShape(Xr_direct(:),p);
      disp(abs(max(Xr_shc(:) - Xr_direct(:))));
      pause(1);
    end
  end

%     for lambda = linspace(0,2*pi,p)
        
%         for m=1:p-1
%           D{m} = [cos(m * lambda) sin(m * lambda);...
%                   -sin(m * lambda) cos(m * lambda)];
%         end
%         vRot = eye(2*p);
%         vRot(2:2*p-1, 2:2*p-1) = blkdiag(D{:});
%         vRot(2*p, 2*p) = cos(p*lambda);
