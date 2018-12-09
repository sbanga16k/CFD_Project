 % Function for solving the intermediate x, y velocity components using the 
% momentum equation employing primite variable formulation
function [pseudo_u, pseudo_v] = pseudo_vel_calc(u, v, ...
    u_bot_nozzles, u_top_nozzles, delta_x, delta_y, delta_t, mu, rho, g_x, g_y)
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
            
            % Solving for u at all nodes except ghost nodes and left wall nodes
            if i < size(u,2)
                % Using FTFS with 1st order upwind for u*du/dx
                if u(j, i) < 0
                    advection_terms = u(j, i) * (u(j, i) - u(j, i-1))/delta_x;
                else
                    advection_terms = u(j, i) * (u(j, i+1) - u(j, i))/delta_x;
                end
                
                % Using FTFS with 1st order upwind for v*d/dy
                if v(j, i) < 0
                    advection_terms = advection_terms + ...
                        v(j, i) * (u(j,i) - u(j-1,i))/(delta_y);
                else
                    advection_terms = advection_terms + ...
                        v(j, i) * (u(j+1,i) - u(j,i))/(delta_y);
                end
                
%                 % TEST
%                 du2_dx = (0.5*0.5*(u(j,i + 1) + u(j,i))*(u(j,i + 1) + u(j,i))... 
%                     - 0.5*0.5*(u(j,i) + u(j,i - 1))*(u(j,i) + u(j,i - 1)))/delta_x;
%                 duv_dx = (0.5*0.5*(u(j + 1,i) + u(j,i))*(v(j,i + 1) + v(j,i))...
%                     -0.5*0.5*(u(j,i) + u(j - 1,i))*(v(j - 1,i + 1) + v(j - 1,i)))/delta_x;
%                 advection_terms = (du2_dx + duv_dx);
                d2u_dx2 = (u(j,i+1) - 2*u(j,i) + u(j,i-1))/(delta_x*delta_x);
                d2u_dy2 = (u(j + 1,i) - 2*u(j,i) + u(j - 1,i))/(delta_y*delta_y);
                pseudo_u(j,i) = u(j,i) + delta_t*(...
                    - advection_terms...
                    + (mu/(rho)) * (d2u_dx2 + d2u_dy2)...
                    ) + g_x * 9.81 * delta_t;
            end
            
            % Solving for v at all nodes except ghost nodes and top & bottom wall nodes
            if j < size(v,1)
                % Using FTFS with 1st order upwind for u*dv/dx
                if u(j, i) < 0
                    advection_terms = u(j, i) * (v(j,i) - v(j,i-1))/(delta_x);
                else
                    advection_terms = u(j, i) * (v(j,i+1) - v(j,i))/(delta_x);
                end
                
                % Using FTFS with 1st order upwind for v*dv/dy
                if v(j, i) < 0
                    advection_terms = advection_terms + ...
                        v(j, i) * (v(j, i) - v(j-1, i))/delta_y;
                else
                    advection_terms = advection_terms + ...
                        v(j, i) * (v(j+1, i) - v(j, i))/delta_y;
                end
%                % TEST
%                 dv2_dx = (0.5*0.5*(u(j,i + 1) + u(j,i))*(u(j,i + 1) + u(j,i))... 
%                     - 0.5*0.5*(u(j,i) + u(j,i - 1))*(u(j,i) + u(j,i - 1)))/delta_x;
%                 dvu_dx = (0.5*0.5*(u(j + 1,i) + u(j,i))*(v(j,i + 1) + v(j,i))...
%                     -0.5*0.5*(u(j,i) + u(j - 1,i))*(v(j - 1,i + 1) + v(j - 1,i)))/delta_x;
%                 advection_terms = (du2_dx + duv_dx);
                d2v_dx2 = (v(j,i + 1) - 2*v(j,i) + v(j,i - 1))/(delta_x*delta_x);
                d2v_dy2 = (v(j+1,i) - 2*v(j,i) + v(j-1,i))/(delta_y*delta_y);
                pseudo_v(j,i) = v(j,i) + delta_t*(...
                    - advection_terms...
                    + (mu/(rho*delta_y^2)) * (d2v_dx2 + d2v_dy2)...
                    ) + g_y * 9.81 * delta_t;
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