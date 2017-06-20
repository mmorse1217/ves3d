function hOut = plotb(Xin, C, varargin)
% PLOTB(X,C) - plots the three dimensional surface with shading C.
%
% The points coordinate vector X is assumed to be returned by the function
% BOUNDARY.
%
% SEE ALSO: BOUNDARY, PARDOMAIN

  if ( isa(Xin,'d3Vec') )
    Xin = vecForm(Xin);
  end
  h = [];

  plotUpRate = 3;
  m = (round(sqrt(2*size(Xin,1)/3+1))-1)/2;
  wasHeld = false;
  if ( ~ishold ), cla; else, wasHeld = true; end
    
  nv = size(Xin,2);
  %- finding the options for this function, cleaning in the input arguments
  %to be passed to surf
  localOpts.interp = any(strcmp(varargin,'interp'));
  localOpts.plain  = any(strcmp(varargin,'plain'));
  ind = find(~(strcmp(varargin,'interp') | strcmp(varargin,'plain')));
  varargin = varargin(ind);
   
  %- Computation
  Xin = reshape(Xin,[], 3 * nv);  
  if(localOpts.interp)
    Xin = interpsh(Xin, plotUpRate*m);
    if (nargin>1 && ~isempty(C)), C = interpsh(C, plotUpRate*m);end 
    m = plotUpRate*m;
  end

  hold on;
  for ii=1:nv
    X = Xin(:,3*(ii-1)+1:3*ii);
    
    %- Putting the caps
    [u, v] = parDomain(m);
    u = reshape(u, m+1,2*m);
    v = reshape(v, m+1,2*m);
    u = [pi*ones(1,2*m); u; 0*ones(1,2*m)];
    v = [v; v(1:2,:)];
    
    for jj=1:3
      XCap(:,jj) = sumBasis( shAna(X(:,jj)), 'Ynm', true, u(:), v(:));
    end
    x = reshape(XCap(:,1), m+3, 2*m);
    y = reshape(XCap(:,2), m+3, 2*m);
    z = reshape(XCap(:,3), m+3, 2*m);

    if ( nargin>1 && ~isempty(C) )
      CC = reshape(C(:,ii), m+1, 2*m);
      CC = [mean(CC(1,:)) * ones(1,2*m); CC; mean(CC(end,:)) * ones(1,2*m)];

      h(ii) = surf([x x(:,1)], [y y(:,1)], [z z(:,1)],[CC CC(:,1)], varargin{:});
      shading interp
    else
      h(ii) = surf([x x(:,1)], [y y(:,1)], [z z(:,1)], varargin{:});
    end

    if (nargin == 1 || isempty(C))
      set(h,'FaceLighting','phong',...
            'FaceColor',[120/255 0 0],...
            'EdgeColor',[0 0 0],...
            'FaceAlpha',.85,...
            'BackFaceLighting','lit',...
            'AmbientStrength',.3,'DiffuseStrength',.8,...
            'SpecularStrength',.9,'SpecularExponent',25);
    end
  end

  if (~localOpts.plain)
    Xin = reshape(Xin,[],3*nv)';
    Xin = reshape(Xin,3,[])';
    mins = min(Xin);
    maxs = max(Xin);
    range = maxs - mins;
    cent = mean(Xin) + [-.4*range(1)/2 -max(range) range(3)/4];
    set(gca,'cameraposition',cent);
    light('Position',cent + [2*range(1)/10 0 range(3)]);
    light('Position',cent + [  range(1)/10 0 -range(3)/2]);
    light('Position',cent + [3*range(1)/10 0 range(3)]);
    camlight('headlight'); lighting gouraud;
  end
  material shiny

  %axis vis3d;
  axis tight;
  axis equal;
  axis off;

  caller = dbstack;
  if(length(caller)<1) %commandline call
    drawnow;
  end
  if(~wasHeld), hold off; end
  if(nargout>0), hOut = h; end
