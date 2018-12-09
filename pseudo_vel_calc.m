% Function for solving the intermediate x, y velocity components using the 
% momentum equation employing primite variable formulation
function [pseudo_u, pseudo_v] = pseudo_vel_calc(u, v, ...
    u_bot_nozzles, u_top_nozzles, delta_x, delta_y, delta_t, mu, rho, g_x, g_y)
    % Shape of u: (n+2) x (m+1)
    % Shape of v: (n+1) x (m+2)
    
    % For the ith element in jth row starting from the bottom, 
    % indices given by - u(i+1/2, j) - u(i-1/2, j), v(i, j+1/2) - v(i, j-1/2)

    % Initializing updated intermediate velocity components to zeroes
    pseudo_u = zeros(size(u));
    pseudo_v = zeros(size(v));

    % Iterating over nodes excluding ghost nodes & boundary nodes
    for j=2:size(v,1)     % Rows: 2 to n+1  (y nodes)
        for i = 2:size(u,2)   %Columns: 2 to m+1  (x nodes)
            
            % Solving for u at all nodes except ghost nodes and left wall nodes
            if i < size(u,2)
                pseudo_u(j,i) = u(j,i) + delta_t*(...
                    - 2 * u(j,i) * (u(j,i) - u(j,i-1))/delta_x...
                    - ( u(j, i) * (v(j,i+1)) - v(j,i)...
                    - u(j, i-1) * (v(j,i)) - v(j,i-1) )/(2*delta_x)...
                    + (mu/(rho*delta_x^2)) * (u(j,i+1) - 2*u(j,i) + u(j,i-1))...
                    ) + g_x * delta_t;
            end
            
            % Solving for v at all nodes except ghost nodes and top & bottom wall nodes
            if j < size(v,1)
                pseudo_v(j,i) = v(j,i) + delta_t*(...
                    - 2 * v(j,i) * (v(j,i) - v(j-1,i))/delta_y...
                    - ( v(j, i) * (u(j+1,i)) - u(j,i)...
                    - v(j-1, i) * (u(j,i)) - u(j-1,i) )/(2*delta_y)...
                    + (mu/(rho*delta_y^2)) * (v(j+1,i) - 2*v(j,i) + v(j-1,i))...
                    ) + g_y * delta_t;
            end
            
        end
    end
    
    % Right boundary nodes update for u velocity?? (Currently flux assumed zero)
    pseudo_u(:,end) = pseudo_u(:,end-1);

    % % Ghost nodes update
    
    % Left wall Ghost nodes update for v velocity
    % to have the average to be equal to the value of v at the wall (0 for our case)
    pseudo_v(:,1) = -pseudo_v(:,2);
    
    % Right wall Ghost nodes update for v velocity (assuming zero flux)
    pseudo_v(:,end) = pseudo_v(:,end-1);
    
    % Bottom wall Ghost nodes update for u velocity
    pseudo_u(1,:) = 2*u_bot_nozzles - pseudo_u(2,:);
    
    % Top wall Ghost nodes update for u velocity
    pseudo_u(end,:) = 2*u_top_nozzles - pseudo_u(end-1,:);
    
end