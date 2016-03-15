% Use U ~1 and theta ~90 degrees and a method

function plot_corner_flow(U, theta_0, method)

    [corner_flow_func] = AP_make_corner_flow(theta_0, U, 'method', method);

    grid = 0.1:0.5:10;
    grid_len = length(grid);

    [x,y] = meshgrid(grid,grid);
    u_x = zeros(grid_len,grid_len);
    u_y = zeros(grid_len,grid_len);
    for i = 1:grid_len
        for j = 1:grid_len
            [V] = corner_flow_func([x(i,j); y(i,j); 3], 10);
            u_x(i, j) = V(1);
            u_y(i,j) = V(2);
        end
    end
    
    quiver(x, y, u_x, u_y)

end

