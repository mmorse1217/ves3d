function PlotShape(X,p,C)

  if(nargin<3), C = []; end
  
  cla;
  hold on;
  for ii = 1:size(X,2)
    Xi = reshape(X(:,ii),[],3);
    
    x = reshape(Xi(:,1), 2*p, p+1)'; 
    y = reshape(Xi(:,2), 2*p, p+1)'; 
    z = reshape(Xi(:,3), 2*p, p+1)'; 
  
    if ( ~isempty(C) )
      cc = reshape(abs(C(:,ii)), 2*p,p+1)'/max(abs(C(:,ii)));
      surf([x x(:,1)], [y y(:,1)], [z z(:,1)], [cc cc(:,1)]);  
      colormap(summer);
    else
      surf([x x(:,1)], [y y(:,1)], [z z(:,1)]);  
      colormap([1 0.3 0.3]);
    end
    
    shading faceted;
    %camlight; lighting gouraud; lightangle(100,0)
    %set(gca,'CameraPosition',[10 -10 3]);
    axis equal;
  end

  hold off;
  axis off;