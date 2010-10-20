function [Rall Rcell]= GenerateShRotMats(p)

  addpath ../../Ves3DMat/src/
  addpath ../../Ves3DMat/util/

  Rall = [];
  Rcell = cell(p+1,1);
  R = cell(p+1,1);
  
  u = parDomain(p);
  for trgIdx = 1:p+1 
    theta = u(trgIdx)-pi;
    
    for n=0:p-1
      R{n+1} = zeros(2*n+1);
    end
    R{p+1} = zeros(2*p);

    for m=1:2*p
      shcIn = zeros(p+1,2*p);
      shcIn(floor(m/2)+1:p+1,m) =1;
      f = shSynReal(shcIn(:));

      fRot = movePole(f,theta,0);
      shcOut = reshape(shAnaReal(fRot(:)),p+1,2*p);
      
      for n=floor(m/2):p
        lim = 2*n +1;
        if(n==p), lim = 2*p;end
        R{n+1}(:,m) = shcOut(n+1,1:lim)';
      end
    end
    % for n=0:p
%       lim = 2*n+1;
%       if(n==p), lim = 2*n; end;
%       for m=1:lim        
%         shcIn = zeros(p+1,2*p);     
%         shcIn(n+1,m) = 1;             
%         f = shSynReal(shcIn(:));
        
%         fRot = movePole(f,theta,0);
%         shcOut = reshape(shAnaReal(fRot(:)),p+1,2*p);
        
%         R{n+1}(:,m) = shcOut(n+1,1:lim)';
%       end
%     end
    
    if ( nargout == 0 )
      X = boundary(p,'dumbbell');
      X = reshape(X.cart.vecForm(),[],3);
      shc = shAnaReal(X);
      for ii=1:3
        temp = reshape(shc(:,ii),p+1,[]);
        shcTemp = zeros(p+1,2*p);
        for jj=0:p
          lim = 2*jj+1;
          if(jj==p), lim = 2*jj;end
          shcTemp(jj+1,1:lim)= R{jj+1} * temp(jj+1,1:lim)';
        end
        shcRot(:,ii) = shcTemp(:);
      end
      X = shSynReal(shcRot);
      plotb(X);
      pause(.1);
    end
    
    for n=0:p  
      Rall = [Rall; R{n+1}(:)];
    end
    Rcell{trgIdx} = R;
  end