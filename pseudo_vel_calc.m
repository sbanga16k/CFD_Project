 % Function for solving the intermediate x, y velocity components using the 
% momentum equation employing primite variable formulation
function [pseudo_u, pseudo_v] = pseudo_vel_calc(u, v, inletLocations,recirculation_uVel,...
    u_bot_nozzles, u_top_nozzles, dx, dy, dt, mu, density, g_x, g_y,inletPresent)
    % Shape of u: (n+2) x (m+1)
    % Shape of v: (n+1) x (m+2)
    
    % For the ith element in jth row starting from the bottom, 
    % indices given by - u(i+1/2, j) - u(i-1/2, j), v(i, j+1/2) - v(i, j-1/2)

    % Initializing updated intermediate velocity components to zeroes
    pseudo_u = u;
    pseudo_v = v;

%     % ORIGINAL
%     % Iterating over nodes excluding ghost nodes & boundary nodes
%     for j=2:size(v,1)     % Rows: 2 to n+1  (y nodes)
%         for i = 2:size(u,2)   %Columns: 2 to m+1  (x nodes)
%             u_curr = u(j,i); 
%             v_curr = v(j,i);
%             
%             u_top = u(j+1,i);
%             u_left = u(j,i-1);
%             u_bot = u(j-1,i);
%             if i < size(u,2)
%                 u_right = u(j,i+1);
%             end
%             
%             if j < size(v,1)
%                 v_top = v(j+1,i);
%             end
%             v_left = v(j,i-1);
%             v_bot = v(j-1,i);
%             v_right = v(j,i+1);
%             
%             if i == size(u,2)-1
%                 u_right = u_curr;   % Zero normal gradient BC
%             end
% 
%             if j == 2               % No slip BC
%                 u_bot = -u_curr;
%             elseif j == size(v,1)   % No slip BC
%                 u_top = -u_curr;
%             end
%             
%             if i == 2
%                 v_left = -v_curr;   % Inlet BC
%             elseif i == size(u,2)
%                 v_right = v_curr;   % Zero grad BC
%             end
%                 
%             % Solving for u at all nodes except ghost nodes and left wall nodes
%             if i < size(u,2)
%                 % u*du/dx
%                 advection_terms = u_curr * (u_right - u_left)/(2*dx);
%                 
%                 % Original eqn, v*du/dy
%                 v_at_u_node = (v_bot + v_curr + v(j-1,i+1) + v_right)/4;     % v-velocity at the location of u(j,i)
%                 advection_terms = advection_terms + ...
%                     v_at_u_node*(u_top - u_bot)/(2*dy);
% 
%                 d2u_dx2 = (u_right - 2*u_curr + u_left)/(dx*dx);
%                 d2u_dy2 = (u_top - 2*u_curr + u_bot)/(dy*dy);
% 
%                 pseudo_u(j,i) = u_curr + dt*( -advection_terms...
%                     + (mu/density) * (d2u_dx2 + d2u_dy2) ) + g_x * dt;
%             end
%             
%             % Solving for v at all nodes except ghost nodes and top & bottom wall nodes
%             if j < size(v,1)
%                 
%                 % Original eqn, u*dv/dx
%                 u_at_v_node = (u_left + u_curr + u(j+1,i-1) + u_top)/4;      % u-velocity at the location of v(j,i)
%                 
%                 advection_terms = u_at_v_node*(v_right - v_left)/(2*dx);
%                 
%                 % v*dv/dy
%                 advection_terms = advection_terms + ...
%                     v_curr*(v_top - v_bot)/(2*dy);
%                 
%                 d2v_dx2 = (v_right - 2*v_curr + v_left)/(dx*dx);
%                 d2v_dy2 = (v_top - 2*v_curr + v_bot)/(dy*dy);
%                 
%                 pseudo_v(j,i) = v_curr + dt*( -advection_terms...
%                     + (mu/density) * (d2v_dx2 + d2v_dy2) ) + g_y * dt;
%             end
%             
%         end
%     end
%     
%     % Right boundary nodes update for u velocity (Zero normal gradient BC)
%     pseudo_u(:,end) = pseudo_u(:,end-1);
% 
%     % % Ghost nodes update
%     
%     % Left wall Ghost nodes update for v velocity
%     % to have the average to be equal to the value of v at the wall (Inlet BC)
%     pseudo_v(:,1) = -pseudo_v(:,2);
%     
%     % Right wall Ghost nodes update for v velocity (Zero normal gradient BC)
%     pseudo_v(:,end) = pseudo_v(:,end-1);
%     
%     % TEMP - CHECKING WITHOUT RECIRCULATION INLETS
% %     % Bottom wall Ghost nodes update for u velocity
% %     pseudo_u(1,:) = 2*u_bot_nozzles - pseudo_u(2,:);
% %     
% %     % Top wall Ghost nodes update for u velocity
% %     pseudo_u(end,:) = 2*u_top_nozzles - pseudo_u(end-1,:);
% 
%     % Bottom wall Ghost nodes update for u velocity
%     pseudo_u(1,:) = -pseudo_u(2,:);
%     
%     % Top wall Ghost nodes update for u velocity
%     pseudo_u(end,:) = -pseudo_u(end-1,:);
    
    
    
    
    
    
    
    % NEW
    numRecirculationInlets = size(inletLocations,1);
    % x-momentum equation
    xPos = 0;
    numUNodes_x = size(u,2);
    numUNodes_y = size(u,1);
    for j = 2:numUNodes_y - 1
        inletLocIter = 1;
        for i = 2:numUNodes_x - 1
            xPos = (i - 1)*dx;                                              % x-coordinate of u-velocity node
            if inletLocIter < numRecirculationInlets && xPos > inletLocations(inletLocIter, 2)
                inletLocIter = inletLocIter + 1;
            end

            % u-velocities at the walls of the CV
%             u_rightWall = u(j,i);
%             u_leftWall = u(j,i - 1);
%             u_topWall = 0.25*(u_leftWall + u(j + 1, i -1) + u(j + 1,i) + u_rightWall);
%             u_bottomWall = 0.25*(u(j - 1, i -1) + u_leftWall + u_rightWall + u(j - 1,i));
            % Neighboring nodes of u(j,i)
            u_curr = u(j,i);
            u_right = u(j,i + 1);
            u_left = u(j,i - 1);
            u_top = u(j + 1,i);
            u_bottom = u(j - 1,i);
            % Checking for BCs (TODO - inlet)
            % BOTTOM (no slip)
            if j == 2
                u_bottom = -u_curr;
                % Checking for an inlet
                if inletPresent && xPos > inletLocations(inletLocIter,1) && xPos > inletLocations(inletLocIter,2)
                    u_bottom = 2*recirculation_uVel - u_curr;
                end
            end
            % TOP (no slip)
            if j == numUNodes_y - 1
                u_top = -u_curr;
                % Checking for an inlet
                if inletPresent && xPos > inletLocations(inletLocIter,1) && xPos > inletLocations(inletLocIter,2)
                    u_top = 2*recirculation_uVel - u_curr;
                end
            end
            % RIGHT (zero gradient)
            if i == numUNodes_x - 1
                u_right = u_curr;
            end
            u_rightWall = 0.5*(u_curr + u_right);
            u_leftWall = 0.5*(u_curr + u_left);
            u_topWall = 0.5*(u_curr + u_top);
            u_bottomWall = 0.5*(u_curr + u_bottom);

            % u-velocities at the neighboring pressure nodes
%             u_curr_pNode = u(j,i);
%             u_left_pNode = u_leftWall;
%             u_right_pNode = u(j, i + 1);
%             u_top_pNode = u(j + 1,i);
%             u_bottom_pNode = u(j - 1,i);

%             u_curr_pNode = u(j,i);
%             u_left_pNode = u(j,i - 1);
%             u_right_pNode = u(j, i + 1);
%             u_top_pNode = u(j + 1,i);
%             u_bottom_pNode = u(j - 1,i);

            % Velocity gradients (Advective terms)
            du_dx = (u_rightWall - u_leftWall)/dx;
            du_dy = (u_topWall - u_bottomWall)/dy;
            
            % Viscous terms
            d2u_dx2 = (u_right - 2*u_curr + u_left)/(dx*dx);
            d2u_dy2 = (u_top - 2*u_curr + u_bottom)/(dy*dy);
            
%             v_curr_pNode = v(j,i);
            v_curr_pNode = 0.25*(v(j,i) + v(j,i + 1) + v(j - 1,i + 1) + v(j - 1,i));
            advectionTerm = -u_curr*du_dx - v_curr_pNode*du_dy;
            viscousTerm = (mu/density)*(d2u_dx2 + d2u_dy2);
            
            pseudo_u(j,i) = u_curr + dt*(advectionTerm + viscousTerm);
        end
    end
    % Updating the u-velocity ghost nodes (TODO - recirculation inlets)
    % Bottom wall, no slip
    pseudo_u(1,:) = -pseudo_u(2,:);
    % Top wall, no slip
    pseudo_u(end,:) = -pseudo_u(end - 1,:);
    % Right wall (zero normal gradient)
    pseudo_u(:,end) = pseudo_u(:,end - 1);
    
    % y-momentum equation
    numVNodes_x = size(v,2);
    numVNodes_y = size(v,1);
    for j = 2:numVNodes_y - 1
        for i = 2:numVNodes_x - 1
            % v-velocities at the walls of the CV
%             v_topWall = v(j,i);
%             v_bottomWall = v(j - 1,i);
%             v_rightWall = 0.25*(v_bottomWall + v_topWall + v(j,i + 1) + v(j - 1,i + 1));
%             v_leftWall = 0.25*(v(j - 1,i - 1)  + v(j,i - 1) + v_topWall + v_bottomWall);

            % Neighboring nodes of v(j,i)
            v_curr = v(j,i);
            v_top = v(j + 1,i);
            v_bottom = v(j - 1,i);
            v_left = v(j,i - 1);
            v_right = v(j,i + 1);
            
            v_topWall = 0.5*(v_curr + v_top);
            v_bottomWall = 0.5*(v_curr + v_bottom);
            v_rightWall = 0.5*(v_curr + v_right);
            v_leftWall = 0.5*(v_curr + v_left);

            % u-velocities at the neighboring pressure nodes
%             v_curr_pNode = v(j,i);
%             v_left_pNode = v(j,i - 1);
%             v_right_pNode = v(j,i + 1);
%             v_top_pNode = v(j + 1,i);
%             v_bottom_pNode = v(j - 1,j);

            % Checking for BCs
            % LEFT (inlet)
            if i == 2
                v_left = -v_curr;
            end
            % RIGHT (zero normal gradient)
            if i == numVNodes_x - 1
                v_right = v_curr;
            end

%             v_curr_pNode = v(j,i);
%             v_left_pNode = v(j,i - 1);
%             v_right_pNode = v(j,i + 1);
%             v_top_pNode = v(j + 1,i);
%             v_bottom_pNode = v(j - 1,j);
            
            % Velocity gradients (Advective terms)
            dv_dx = (v_rightWall - v_leftWall)/dx;
            dv_dy = (v_topWall - v_bottomWall)/dy;
            
            % Viscous terms
            d2v_dx2 = (v_right - 2*v_curr + v_left)/(dx*dx);
            d2v_dy2 = (v_top - 2*v_curr + v_bottom)/(dy*dy);
            
%             u_curr_pNode = u(j,i);
            u_curr_pNode = 0.25*(u(j,i) + u(j,i - 1) + u(j + 1,i - 1) + u(j + 1,i));
            advectionTerm = -u_curr_pNode*dv_dx - v_curr*dv_dy;
            viscousTerm = (mu/density)*(d2v_dx2 + d2v_dy2);
            
            pseudo_v(j,i) = v_curr + dt*(advectionTerm + viscousTerm + g_y);
        end
    end
    % Updating the v-velocity ghost nodes (TODO - recirculation inlets)
    % Left wall (Average is zero)
    pseudo_v(:,1) = -pseudo_v(:,2);
    % Right wall (zero normal gradient)
    pseudo_v(:,end) = pseudo_v(:,end - 1);

end