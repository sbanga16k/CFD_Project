% Function for updating x, y velocity components for next time step using the 
% momentum equation employing primite variable formulation
function [u_new, v_new] = update_vel(pseudo_u, pseudo_v, p, delta_t, dx,dy,...
    u_bot_nozzles, u_top_nozzles, density)
    % Shape of u: (n+2) x (m+1)
    % Shape of v: (n+1) x (m+2)
    
    % For the ith element in jth row starting from the bottom, 
    % indices given by - u(i+1/2, j) - u(i-1/2, j), v(i, j+1/2) - v(i, j-1/2)
    
    % Initializing updated intermediate velocity components to zeroes
    u_new = pseudo_u;
    v_new = pseudo_v;
    
    % Calculate correct velocity components using pressure at next time
    % step for all nodes except boundary & ghost nodes
    for j=2:size(pseudo_v,1)         % Rows: 2 to n+1  (y nodes)
        for i = 2:size(pseudo_u,2)   % Columns: 2 to m+1  (x nodes)
%             delta_t*(p(j,i) - p(j,i-1))/(density*9.81)
            u_new(j,i) = pseudo_u(j,i) - delta_t*(p(j,i) - p(j,i-1))/(density*dx);
            v_new(j,i) = pseudo_v(j,i) - delta_t*(p(j,i) - p(j-1,i))/(density*dy);
            
        end
    end

    % Right boundary nodes update for u velocity?? (Currently flux assumed zero)
    u_new(:,end) = u_new(:,end-1);

    % % Ghost nodes update
    
    % Left wall Ghost nodes update for v velocity
    % to have the average to be equal to the value of v at the wall (0 for our case)
    v_new(:,1) = -v_new(:,2);
    
    % Right wall Ghost nodes update for v velocity (assuming zero flux)
    v_new(:,end) = v_new(:,end-1);
    
    % Bottom wall Ghost nodes update for u velocity
    u_new(1,:) = 2*u_bot_nozzles - u_new(2,:);
    
    % Top wall Ghost nodes update for u velocity
    u_new(end,:) = 2*u_top_nozzles - u_new(end-1,:);
    
end