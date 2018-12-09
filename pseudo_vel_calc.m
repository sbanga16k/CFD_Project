 % Function for solving the intermediate x, y velocity components using the 
% momentum equation employing primite variable formulation
function [pseudo_u, pseudo_v] = pseudo_vel_calc(u, v, ...
    u_bot_nozzles, u_top_nozzles, dx, dy, dt, mu, density, g_x, g_y)
    % Shape of u: (n+2) x (m+1)
    % Shape of v: (n+1) x (m+2)
    
    % For the ith element in jth row starting from the bottom, 
    % indices given by - u(i+1/2, j) - u(i-1/2, j), v(i, j+1/2) - v(i, j-1/2)

    % Initializing updated intermediate velocity components to zeroes
    pseudo_u = u;
    pseudo_v = v;

    % Iterating over nodes excluding ghost nodes & boundary nodes
    for j=2:size(v,1)     % Rows: 2 to n+1  (y nodes)
        for i = 2:size(u,2)   %Columns: 2 to m+1  (x nodes)
            u_curr = u(j,i); 
            v_curr = v(j,i);
            
            u_top = u(j+1,i);
            u_left = u(j,i-1);
            u_bot = u(j-1,i);
            if i < size(u,2)
                u_right = u(j,i+1);
            end
            
            if j < size(v,1)
                v_top = v(j+1,i);
            end
            v_left = v(j,i-1);
            v_bot = v(j-1,i);
            v_right = v(j,i+1);
            
            if i == size(u,2)-1
                u_right = u_curr;   % Zero normal gradient BC
            end

            if j == 2               % No slip BC
                u_bot = -u_curr;
            elseif j == size(v,1)   % No slip BC
                u_top = -u_curr;
            end
            
            if i == 2
                v_left = -v_curr;   % Inlet BC
            elseif i == size(u,2)
                v_right = v_curr;   % Zero grad BC
            end
                
            % Solving for u at all nodes except ghost nodes and left wall nodes
            if i < size(u,2)
                % u*du/dx
                advection_terms = u_curr * (u_right - u_left)/(2*dx);
                
                % Original eqn, v*du/dy
                v_curr = (v_bot + v_curr + v(j-1,i+1) + v_right)/4;     % v-velocity at the location of u(j,i)
                advection_terms = advection_terms + ...
                    v_curr*(u_top - u_bot)/(2*dy);

                d2u_dx2 = (u_right - 2*u_curr + u_left)/(dx*dx);
                d2u_dy2 = (u_top - 2*u_curr + u_bot)/(dy*dy);

                pseudo_u(j,i) = u_curr + dt*( -advection_terms...
                    + (mu/density) * (d2u_dx2 + d2u_dy2) ) + g_x * dt;
            end
            
            % Solving for v at all nodes except ghost nodes and top & bottom wall nodes
            if j < size(v,1)
                
                % Original eqn, u*dv/dx
                u_curr = (u_left + u_curr + u(j+1,i-1) + u_top)/4;      % u-velocity at the location of v(j,i)
                
                advection_terms = u_curr*(v_right - v_left)/(2*dx);
                
                % v*dv/dy
                advection_terms = advection_terms + ...
                    v_curr*(v_top - v_bot)/(2*dy);
                
                d2v_dx2 = (v_right - 2*v_curr + v_left)/(dx*dx);
                d2v_dy2 = (v_top - 2*v_curr + v_bot)/(dy*dy);
                
                pseudo_v(j,i) = v_curr + dt*( -advection_terms...
                    + (mu/density) * (d2v_dx2 + d2v_dy2) ) + g_y * dt;
            end
            
        end
    end
    
    % Right boundary nodes update for u velocity (Zero normal gradient BC)
    pseudo_u(:,end) = pseudo_u(:,end-1);

    % % Ghost nodes update
    
    % Left wall Ghost nodes update for v velocity
    % to have the average to be equal to the value of v at the wall (Inlet BC)
    pseudo_v(:,1) = -pseudo_v(:,2);
    
    % Right wall Ghost nodes update for v velocity (Zero normal gradient BC)
    pseudo_v(:,end) = pseudo_v(:,end-1);
    
    % Bottom wall Ghost nodes update for u velocity
    pseudo_u(1,:) = 2*u_bot_nozzles - pseudo_u(2,:);
    
    % Top wall Ghost nodes update for u velocity
    pseudo_u(end,:) = 2*u_top_nozzles - pseudo_u(end-1,:);
    
end