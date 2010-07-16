function err = GenerateBmark_TwoVes_Imp(precision)
  
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
    saveThis(VES3D_PATH, 'x0.txt', ves, precision);
    subplot(1,2,1); PlotShape(ves);
    
    Nv = length(ves); 
    Fb = cell(Nv,1); Ws = cell(Nv,1);
    Fsig = cell(Nv,1); [Fsig{:}] = deal(zeros(6*m*(m+1),1));
    FI = cell(Nv,1); [FI{:}] = deal(zeros(6*m*(m+1),1));
    tension = cell(Nv,1);
    MatVecRes = cell(Nv,1);
    
    for k = 1:Nv
      [H, K, W, nor, GradS, DivS] = Curvatures(ves{k});  
      Ws{k} = W; 
      Fn = -DivS(GradS(H)) - 2*H.*(H.^2 - K);   
      Fb{k} = kappa*repmat(Fn, 3, 1).*nor;                  
    end
    
    for p = 1:Nv
      for q = setdiff(1:Nv, p)
        FI{p} = FI{p} + DirectStokesNew(ves{q}, Fb{q}+Fsig{q}, Ws{q}, ves{p});
      end
    end
    
    for k = 1:Nv
      GS = kernelS(ves{k}, Ws{k}); 
      [H, K, W, nor, GradS, DivS, E, F, G] = Curvatures(ves{k});  
      FI{k} = FI{k} + shear*[ves{k}(1+2*end/3:end); zeros(4*m*(m+1),1)];
      
      MatVecRes{k} = TimeSchemeII(ves{k}, H, K, DivS, GradS, nor, E, F, G, W, GS, ...
                               kappa, ts);
      
      %------ local constraint
      Fs = @(sig) GS*(GradS(sig) + 2*repmat(sig.*H,3,1).*nor);  
      L =  @(sig) DivS(Fs(sig));                                     
      [tension{k},FLAG,RELRES,ITER] = gmres(L, -DivS(FI{k}+GS*Fb{k}), ...
                                            [], 10^-8, 2*m*(m+1));
      
      FI{k} = FI{k} + Fs(tension{k});
      
      %------- time marching
      [Xn,FLAG,RELRES,ITER] = gmres(@TimeSchemeII, ves{k}+ts*FI{k}, ...
                                    [], 1e-8, numel(ves{k}),[],[],[], ...
                                    H, K, DivS, GradS, nor, E, F, G, W, GS, kappa, ts); 
      ves{k} = Xn; 
    end
    
    saveThis(VES3D_PATH, 'matvec.txt', MatVecRes, precision);
    saveThis(VES3D_PATH, 'tension.txt', tension, precision);
    saveThis(VES3D_PATH, 'xnew.txt', ves, precision);
    
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
  
  name = [path filesep  'Bmark_TwoVes_Imp_p' num2str(m) '_' precision '_' name];
  disp(['Saving ' name]);

  save(name,'array','-ascii'); 
%   fid = fopen(name,'w');
%   fwrite(fid, array,precision);
%   fclose(fid);


  
