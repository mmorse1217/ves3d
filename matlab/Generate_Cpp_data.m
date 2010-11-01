function Generate_Cpp_data(p, precision, varargin)
  
  addpath ../../Ves3DMat/src/
  addpath ../../Ves3DMat/util/

  for ii = 1:length(varargin)
    
    switch varargin{ii}
     case 'SpHarmRot'
      R = GenerateShRotMats(p);
      R = R(:);
      fileName = ['SpHarmRotMats_p' num2str(p) '_' precision];
      saveData([fileName '.bin'], R, precision);  
      save([fileName '.txt'],'R','-ascii'); 
      
     case 'LegendreTransform'
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

      fileName = ['legTrans' num2str(p)  '_' precision '.txt'];
      save(fileName,'Lmat','-ascii'); 

      fileName = ['legTransInv' num2str(p)  '_' precision '.txt'];
      save(fileName,'LImat','-ascii'); 

      fileName = ['d1legTrans' num2str(p)  '_' precision '.txt'];
      save(fileName,'Hmat','-ascii'); 

      fileName = ['d2legTrans' num2str(p)  '_' precision '.txt'];
      save(fileName,'Wmat','-ascii'); 


     case 'DumbbellShape'
      S = boundary(p,'dumbbell');
      S.resample(p);
      S.plot;
      x = reshape(S.cart.x,p+1,[])';
      y = reshape(S.cart.y,p+1,[])';
      z = reshape(S.cart.z,p+1,[])';

      X = [x(:); y(:); z(:)];
      fileName = ['dumbbell_' num2str(p) '_' precision '.txt'];
      save(fileName,'X','-ascii'); 

     case 'DirectRotation'
      R = GenerateDirectRotMat(p);
      R = R(:);
      fileName = ['all_rot_mats_' num2str(p)  '_' precision '.txt'];
      save(fileName,'R','-ascii'); 

     case 'SingularIntegralWeights'
      sqw = SingularWeights(p);
      sqw = reshape(sqw,p+1,[]);
      sqw = flipud(sqw);
      sqw = sqw';
      sqw = sqw(:);
      fileName = ['sing_quad_weights_' num2str(p)  '_' precision '.txt'];
      save(fileName,'sqw','-ascii'); 

     case 'IntegrationWeights'
      [trash qw] = grule(p+1);
      qw = repmat(qw(:)',2*p,1);
      qw = qw(:)*pi/p;
      u = parDomain(p);
      u = reshape(u,p+1,[])';
      qw = qw./sin(u(:));
      fileName = ['quad_weights_' num2str(p)  '_' precision '.txt'];
      save(fileName,'qw','-ascii'); 

     case 'WSpherical'
      u = parDomain(p);
      u = reshape(u,p+1,[])';
      w_sph = sin(u(:));

      fileName = ['w_sph_' num2str(p)  '_' precision '.txt'];
      save(fileName,'w_sph','-ascii'); 

    end      
  end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shape = 'two_ellipse_';
% 
% [u,v] = parDomain(p); 
% a = 2;
% X =  sin(u).*cos(v);     
% Y = sin(u).*sin(v);
% Z = a*cos(u);
% 
% cart(1) = d3Vec([X+2;Y;Z-.7]);
% cart(2) = d3Vec([X-2;Y;Z+.7]);
% 
% for ii=1:length(cart)
%   cart(ii).x = reshape(cart(ii).x,p+1,[])';
%   cart(ii).x = cart(ii).x(:);
%   cart(ii).y = reshape(cart(ii).y,p+1,[])';
%   cart(ii).y = cart(ii).y(:);
%   cart(ii).z = reshape(cart(ii).z,p+1,[])';
%   cart(ii).z = cart(ii).z(:);
% end
% 
% X = [cart(1).vecForm; cart(2).vecForm];
% fileName = [shape num2str(p)];
% save(fileName,'X','-ascii'); 

    %%% 

% NV = 600;
% zMax = 80;
% rMax = 20;
% minDist = 5;
% 
% [x,y,z] = pol2cart(0,0,zMax/2);
% centers = [x y z];
% while size(centers,1)<NV
%   z = rand*zMax;
%   r = rand*rMax;
%   t = rand*2*pi;
%   
%   [x,y,z] = pol2cart(t,r,z);
%   d = centers - repmat([x y z],size(centers,1),1);
%   d = sqrt(min(dot(d,d,2)));
%   if(d>minDist)
%     centers = [centers;[x y z]];
%   end
%   
%   disp(size(centers,1));
% end
% centers = [centers(:,3) centers(:,2) centers(:,1)];
% scatter3(centers(:,1),centers(:,2),centers(:,3));
% 
% centers = centers';
% centers = centers(:);
% fileName = 'centers_tube_600.txt';
% save(fileName,'centers','-ascii'); 

% [r t] = meshgrid(0:10:20,0:pi/4:(2*pi-pi/4));
% r = r(:,2:end); r=r(:);
% t = t(:,2:end); t=t(:);
% r = [0;r];
% t = [0;t];
% 
% [x y] = pol2cart(t,r);
% plot(x,y,'o');
% c = [zeros(size(x)) y x];
% 
% fileName = 'pepperoni.txt';
% save(fileName,'c','-ascii'); 

% NV = 10;
% side = 8;
% minDist = 2;

% centers = [0 0 0];
% cc = S.cart;
% while size(centers,1)<NV
%   xx = rand(1,3)*side;
%   xx(1) = xx(1)/side;
    
%   d = centers - repmat(xx,size(centers,1),1);
%   d = sqrt(min(dot(d,d,2)));
%   if(d>minDist)
%     centers = [centers;xx];
%     cc(end+1).x = S.cart.x + xx(1);
%     cc(end).y = S.cart.y + xx(2);
%     cc(end).z = S.cart.z + xx(3);    
%   end
    
%   disp(size(centers,1));
% end
% PlotShape(cc,12)

