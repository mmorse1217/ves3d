clear all;clc;clf

p = 6;
nv = 2;
np = 2 * p * (p + 1);
nSnapShot = 100;

xname    = './Surfaces.out';
gridname = './EvaluationPoints.out';
velname  = '../test/VelocityField.out';

%- Reading the configuration
if ( exist(xname) )
  xv = ReadAscii(xname);
  xv = reshape(xv,[],nv * nSnapShot);
  xv = RowToColMajor(xv,p);
  xv = reshape(xv,[], nSnapShot);
end

%- Generating the grid points
gridsize = {50,20,20};
gridlim  = [-3 3;-1.5 1.5;-1.5 1.5];
gridSpec = {linspace(gridlim(1,1),gridlim(1,2),gridsize{1}),
            linspace(gridlim(2,1),gridlim(2,2),gridsize{2}),
            linspace(gridlim(3,1),gridlim(3,2),gridsize{3})};           

[x y z] = meshgrid(gridSpec{:});
XX = [x(:);y(:);z(:)];
save(gridname,'XX','-ascii'); 
gridsize = mat2cell(size(x),1,[1 1 1]);

%- Reading streamline data
if ( exist(velname) )
  XAll = ReadAscii(gridname);
  UAll = reshape(ReadAscii(velname),[], nSnapShot);
  
  s_gridsize = {2,6,4};
  s_gridSpec = {linspace(gridlim(1,1),gridlim(1,2),s_gridsize{1}),
                linspace(.8 * gridlim(2,1),.8 * gridlim(2,2),s_gridsize{2}),
                linspace(.8 * gridlim(3,1),.8 * gridlim(3,2),s_gridsize{3})};           
  
  [x y z] = meshgrid(s_gridSpec{:});
  xs = [squeeze(x(:,end,end/2+1:end)) squeeze(x(:,1,1:end/2))];
  ys = squeeze(y(:,1,:));
  zs = squeeze(z(:,1,:));

  for ii=1:65
    disp(ii);
    XX = reshape(XAll,[],3);
    UU = reshape(UAll(:,ii),[],3);

    xx = reshape(XX(:,1), gridsize{:});
    yy = reshape(XX(:,2), gridsize{:});
    zz = reshape(XX(:,3), gridsize{:});

    uu = reshape(UU(:,1), gridsize{:});
    vv = reshape(UU(:,2), gridsize{:});
    ww = reshape(UU(:,3), gridsize{:});

%   hold on;
%   quiver3(xx,yy,zz,uu,vv,ww,.5);
%   axis(reshape(gridlim',1,[]));
%   axis on;
%   hold off;

    clf;
    hold on;
    Xplot = reshape(xv(:,ii),[],nv);
    for jj=1:nv
      Plotb(Xplot(:,jj),p);
    end
    vertices = stream3(xx,yy,zz,uu,vv,ww, xs, ys, zs);
    streamtube(vertices,.05);
    view(3);
    axis tight;
    colormap([120/255 0 0]);
    shading interp;
    camlight; lighting gouraud;
    lightangle(100,0)
    lightangle(-100,0)
    hold off;
    set(gca,'CameraPosition',[5 -15 3]);
    axis([-3 3 -1.6 1.6 -1.6 1.6]);
    %f = figure('visible','off')
    drawnow;
    saveas(gcf,['stplot' num2str(ii,'%03g')],'png');
    pause(.5);
  end
end
