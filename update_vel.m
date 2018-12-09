% Function for updating x, y velocity components for next time step using the 
% momentum equation employing primite variable formulation
function [u_new, v_new] = update_vel(pseudo_u, pseudo_v, p, delta_t, ...
    u_bot_nozzles, u_top_nozzles)
    % Shape of pseudo_u: m x (n+2)
    % Shape of pseudo_v: (m+2) x n
    
    % For the ith element in jth row starting from the bottom, 
    % indices given by - u(i+1/2, j) - u(i-1/2, j), v(i, j+1/2) - v(i, j-1/2)
    
    % Initializing updated intermediate velocity components to zeroes
    u_new = zeros(size(pseudo_u));
    v_new = zeros(size(pseudo_v));
    
    % Calculate correct velocity components using pressure at next time
    % step for all nodes except boundary & ghost nodes
    for j=1:size(v,2)     % Rows: 1 to n  (y nodes)
        for i = 1:size(u,1)   %Columns: 1 to m  (x nodes)
            % For u, j+1 corresponds to cell i,j
            % For v, i+1 corresponds to cell i,j
            
            if (i >= 2)
                u_new(i,j+1) = pseudo_u(j+1,i) - delta_t*(p(j,i+1) - p(j,i));
            end
            
            if (j >= 2)
                v_new = pseudo_v(j,i+1) - delta_t*(p(j+1,i) - p(j,i));
            end
            
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