function PlotShape(X,p)

  cla;
  hold on;
  for ii = 1:size(X,2)
    Xi = reshape(X(:,ii),[],3);
    
    x = reshape(Xi(:,1), 2*p, p+1)'; 
    y = reshape(Xi(:,2), 2*p, p+1)'; 
    z = reshape(Xi(:,3), 2*p, p+1)'; 
    surf([x x(:,1)], [y y(:,1)], [z z(:,1)]);  
    
    camlight; lighting gouraud; lightangle(100,0)
    set(gca,'CameraPosition',[10 -10 3]);
    axis equal;
    colormap([1 0.3 0.3]) 
  end
  hold off;
  axis off;