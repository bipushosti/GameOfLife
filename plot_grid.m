
function plot_grid(result)
  
  grid_dim = size(result)(2);
  iterations = size(result)(1) / grid_dim;
  
  [x, y] = meshgrid(1:grid_dim, 1:grid_dim); 

  surf(x, y, result(1:grid_dim,:)); 
  view([0 -90]);
  axis ([1 grid_dim 1 grid_dim 0 1]);
  title("Iteration: 1");
  hold on 
  drawnow 


  for i = 2:iterations 
    
    pause(0.25);
    row_start = (i - 1) * grid_dim + 1;
    row_end = i * grid_dim;
   
    clf;
    #surf( x, y, F(x,y,t(i)) ); 
    surf(x,y,result(row_start:row_end,:));
    view([0 -90]);
    axis ([1 grid_dim 1 grid_dim 0 1]);
    title(sprintf('Iteration: %d', i))
    drawnow 
    
    
  end 
endfunction