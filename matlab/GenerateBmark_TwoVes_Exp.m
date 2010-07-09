function err = GenerateBmark_TwoVes_Exp(precision)
  
  VES3D_PATH = pwd;
  MATLAB_PATH = '/home/abtin/school/projects/Ves3DMat/Shravan/';
  
  err = [];
  try
    cd(MATLAB_PATH);
    global m;
    m = 12;
    
    %% Initializing and reading data
    main_startup;
    kappa = .01; shear = 0.1; ts = 1;

    %% boundary 
    ves = boundaryMultiple(m);
    saveThis(VES3D_PATH, 'x0.bin', ves, precision);
    subplot(1,2,1); PlotShape(ves);
    
    Nv = length(ves); GS = cell(Nv,1); 
    Fb = cell(Nv,1); W = cell(Nv,1); nor = cell(Nv,1);
    Fsig = cell(Nv,1); [Fsig{:}] = deal(zeros(6*m*(m+1),1));
    FI = cell(Nv,1); [FI{:}] = deal(zeros(6*m*(m+1),1));
    
    for k = 1:Nv
      [H0, K, W0, nor0, GradS0, DivS0] = Curvatures(ves{k});  
      W{k} = W0; nor{k} = nor0; GradS{k} = GradS0; DivS{k} = DivS0; H{k} = H0;
      
      Fn = -DivS{k}(GradS{k}(H{k})) - 2*H{k}.*(H{k}.^2 - K);   
      Fb{k} = kappa*repmat(Fn, 3, 1).*nor{k};                  
    end
    saveThis(VES3D_PATH, 'Fb.bin', Fb, precision);
    
    for p = 1:Nv
      for q = setdiff(1:Nv, p)
        FI{p} = FI{p} + DirectStokesNew(ves{q}, Fb{q}+Fsig{q}, W{q}, ves{p});
      end
    end
      
    for k = 1:Nv        
      G = kernelS(ves{k}, W{k});       
      SFb{k} = G*Fb{k};
      
      F = shear*[ves{k}(1+2*end/3:end); zeros(4*m*(m+1),1)] + G*Fb{k} + FI{k};
      vel{k} = F;
      
      %------ local constraint
      Fs = @(sig) G*(GradS{k}(sig) + 2*repmat(sig.*H{k},3,1).*nor{k});  
      L =  @(sig) DivS{k}(Fs(sig));                                     
      [sig,FLAG,RELRES,ITER] = gmres(L, -DivS{k}(F), [], 10^-8, 2*m*(m+1));
      F = F + Fs(sig);                                                    
      tension{k} = sig;
      
      %------- time marching
      ves{k} = ves{k} + ts*F;                  
    end
      
    saveThis(VES3D_PATH, 'SFb.bin', SFb, precision);
    saveThis(VES3D_PATH, 'vel.bin', vel, precision);
    saveThis(VES3D_PATH, 'tension.bin', tension, precision);
    saveThis(VES3D_PATH, 'xnew.bin', ves, precision);
    subplot(1,2,2); PlotShape(ves);
  catch err
    disp(' There was an error in invoking the Matlab code:');
    disp(err);
  end
  cd(VES3D_PATH);
  
function saveThis(path, name, array, precision)

  if(iscell(array))
    temp = [];
    for ii=1:length(array)
      temp = [temp; array{ii}(:)];
    end
    array = temp;
  end
  
  global m;
  np = 2 * m * ( m + 1);
  %% saving it as row major
  array = reshape(array,np,[]);
  for  ii = 1:size(array,2)
    array(:,ii) = reshape(reshape(array(:,ii),m+1, 2*m)',[],1);
  end
  array = array(:);     
  
  name = [path filesep  'Bmark_TwoVes_Exp_p' num2str(m) '_' precision '_' name];
  disp(['Saving ' name]);
  fid = fopen(name,'w');
  fwrite(fid, array,precision);
  fclose(fid);
  
