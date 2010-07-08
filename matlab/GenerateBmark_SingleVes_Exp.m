function err = GenerateBmark_singleVes_Exp(precision)

  VES3D_PATH = pwd;
  MATLAB_PATH = '/home/abtin/school/projects/Ves3DMat/Shravan/';
  
  err = [];
  try
    cd(MATLAB_PATH);
    global m;
    m = 12;
    
    %% Generating the shape
    [u,v] = parDomain(m); rho = 1 + real(Ynm(2,0,u,v));  
    X = [rho.*sin(u).*cos(v); rho.*sin(u).*sin(v); rho.*cos(u)];
    saveThis(VES3D_PATH, 'x0.bin', X, precision);
    
    %% Initializing and reading data
    main_startup;
    kappa = .01; shear = 0.1; ts = 1;
    
    %% Differentiating 
    [H, K, W, nor, GradS, DivS] = Curvatures(X);
    G = kernelS(X, W);
    
    %% Bending and background flow
    v_inf = shear*[X(1+2*end/3:end); zeros(4*m*(m+1),1)];               
    Fn = - DivS(GradS(H)) - 2*H.*(H.^2 - K);                             
    Fb = kappa*repmat(Fn, 3, 1).*nor; 
    saveThis(VES3D_PATH, 'Fb.bin', Fb, precision);
    vel = G*Fb; 
    saveThis(VES3D_PATH, 'SFb.bin', vel, precision);
    vel = vel + v_inf;
    saveThis(VES3D_PATH, 'vel.bin', vel, precision);
    
    %% Tension
    Fs = @(sig) G*(GradS(sig) + 2*repmat(sig.*H,3,1).*nor);             
    L =  @(sig) DivS(Fs(sig));                                          
    [sig,FLAG,RELRES,ITER] = gmres(L, -DivS(vel), [], 10^-8, 2*m*(m+1));  
    saveThis(VES3D_PATH, 'tension.bin', sig, precision);
    vel = vel + Fs(sig);                                                    

    %% New position
    Xn = X + ts*vel; 
    saveThis(VES3D_PATH, 'xnew.bin', Xn, precision);
    
    subplot(1,2,1); plotb(X,m); alpha(0.55); axis tight; axis off; 
    subplot(1,2,2); plotb(Xn,m); alpha(0.55); axis tight; axis off; 
    
  catch err
    disp(' There was an error in invoking the Matlab code:');
    disp(err);
  end
  cd(VES3D_PATH);
  
function saveThis(path, name, array, precision)

  global m;
  np = 2 * m * ( m + 1);
  %% saving it as row major
  array = reshape(array,np,[]);
  for  ii = 1:size(array,2)
    array(:,ii) = reshape(reshape(array(:,ii),m+1, 2*m)',[],1);
  end
  array = array(:);     
  
  name = [path filesep  'Bmark_SingleVes_Exp_p' num2str(m) '_' precision '_' name];
  disp(['Saving ' name]);
  fid = fopen(name,'w');
  fwrite(fid, array,precision);
  fclose(fid);
  
