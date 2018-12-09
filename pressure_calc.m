% Point Gauss-Siedel method for solving pressure equation using intermediate 
% values of u,v velocities by formulating it as elliptical PDE with provision for
% Successive over-relaxation (SOR)
function [pressure_new] = pressure_calc(pressure, pseudo_u, pseudo_v,...
    p_inlet, inletLocations, p_bot_nozzles,density, ...
    delta_x, delta_y, delta_t, sor_factor, epsilon)
    % Shape of pseudo_u: (n+2) x (m+1)
    % Shape of pseudo_v: (n+1) x (m+2)
    % Shape of pressure: (n+2) x (m+2)
    
    beta = delta_x/delta_y;

    % Intializing new matrix to store updated values of phi for testing
    % convergence of algorithm by computing difference with phi values at
    % previous time step
%     pressure_new = zeros(size(pressure));
    pressure_new = pressure;

    % Stores max error in the matrix for testing algo convergence
    max_error = 100;

    % Recirculation Inlet info
    numRecirculationInlets = size(inletLocations, 1);
    % Iterating over all the interior nodes to update the phi values till the
    % algo converges
    while max_error > epsilon
        for j = 2:size(pseudo_v, 1)
            xPos = delta_x/2;                        % x-coordinate of pressure node
            inletLocIter = 1;
            for i = 2:size(pseudo_u, 2)
                currNodeCoeff = (2.0 / (delta_x*delta_x)) + (2.0 / (delta_y*delta_y));  % Coefficient of P_j,i
                xPos = delta_x/2 + (i - 2)*delta_x;
                if inletLocIter < numRecirculationInlets && xPos > inletLocations(inletLocIter, 2)
                    inletLocIter = inletLocIter + 1;
                end
%                 pressure_new(j,i) = (1 - sor_factor)* pressure(j,i) + ...
%                     (0.5*sor_factor/(1 + beta*beta))*(pressure_new(j, i-1) + ...
%                     pressure(j, i+1) + ...
%                     (beta*beta)*(pressure_new(j-1, i) + pressure(j+1,i))) - ...
%                     0.5/(1 + beta*beta) * (density*delta_x^2/delta_t) * ...
%                     ((pseudo_u(j,i) - pseudo_u(j,i-1))/delta_x + ...
%                     (pseudo_v(j,i) - pseudo_v(j-1,i))/delta_y);

                % NEW
                % Left wall
                p_curr = pressure(j,i);
                if i == 2
                    p_left = 2*p_inlet;
                    currNodeCoeff = currNodeCoeff + 1/(delta_x*delta_x);
                else
                    p_left = pressure_new(j, i-1);
                end
                % Right wall-
                if i == size(pseudo_u, 2)
                    p_right = 0;
                    currNodeCoeff = currNodeCoeff - 1/(delta_x*delta_x);
                else
                    p_right = pressure(j, i + 1);           % Zero gradient condition
                end
                % Bottom wall (TODO - inlet)
                if j == 2
                    % Inlet exists at the bottom wall
                    if xPos > inletLocations(inletLocIter,1) && xPos < inletLocations(inletLocIter,2)
                        p_bottom = 2*p_inlet;      % At the bottom wall inlet
                        currNodeCoeff = currNodeCoeff + 1/(delta_y*delta_y);
                    else
                        p_bottom = 0;       % Zero gradient condition
                        currNodeCoeff = currNodeCoeff - 1/(delta_y*delta_y);
                    end
                else
                    p_bottom = pressure_new(j - 1, i);
                end
                % Top wall (TODO - inlet)
                if j == size(pseudo_v, 1)
                    % Inlet exists at the top wall
                    if xPos > inletLocations(inletLocIter,1) && xPos < inletLocations(inletLocIter,2)
                        p_top = 2*p_inlet;         % At the top wall inlet
                        currNodeCoeff = currNodeCoeff + 1/(delta_y*delta_y);
                    else
                        p_top = 0;              % Zero gradient condition
                        currNodeCoeff = currNodeCoeff - 1/(delta_y*delta_y);
                    end
                else
                    p_top = pressure(j + 1, i);
                end
                new_p_curr = (p_left + p_right)/(delta_x*delta_x) + (p_top + p_bottom)/(delta_y*delta_y)...
                        - density*((pseudo_u(j,i) - pseudo_u(j,i-1))/delta_x + ...
                (pseudo_v(j,i) - pseudo_v(j-1,i))/delta_y)/(delta_t);
                new_p_curr = new_p_curr/currNodeCoeff;
                % SOR
                pressure_new(j,i) = (1 - sor_factor)*pressure(j,i) + sor_factor*new_p_curr;                
            end
        end
        
%         % % Ghost nodes update for pressure
%        
%         % Left wall Ghost nodes update
%         % to have the average to be equal to the value of pressure at inlet (p_inlet)
%         pressure_new(:,1) = 2*p_inlet - pressure_new(:,2);
% 
%         % Right wall Ghost nodes update (assuming zero flux)
%         pressure_new(:,end) = pressure_new(:,end-1);
% 
%         % Bottom wall Ghost nodes update (assuming zero flux)
%         pressure_new(1,:) = pressure_new(2,:);
%         for k = 1:size(p_bot_nozzles, 1)
%             nozzle_inds = p_bot_nozzles(k,:)./delta_x + 1
%             pressure_new(1,nozzle_inds(1):nozzle_inds(2)) = 2*p_inlet - ...
%                 pressure_new(2,nozzle_inds(1):nozzle_inds(2));
%         end
% 
%         % Top wall Ghost nodes update (assuming zero flux)
%         pressure_new(end,:) = pressure_new(end-1,:);
%         for l = 1:size(p_top_nozzles, 1)
%             nozzle_inds = p_top_nozzles(l,:)./delta_x + 1;
%             pressure_new(end,nozzle_inds(1):nozzle_inds(2)) = 2*p_inlet - ...
%                 pressure_new(end-1,nozzle_inds(1):nozzle_inds(2));
%         end

        % Compute error in values between consecutive time steps and sets
        % pressure to be equal to pressure_new
        max_error = max(max(abs(pressure_new(2:end-1,2:end-1) - pressure(2:end-1,2:end-1))));
        pressure = pressure_new;
    end
    % Updating the pressure ghost nodes
    % Left wall (Average should be equal to inlet pressure)
    pressure_new(:,1) = 2*p_inlet*ones(numel(pressure_new(:,2)),1) - pressure_new(:,2);
    
    % Right wall (Zero gradient condition)
    pressure_new(:, end) = pressure_new(:,end - 1);
    
    % Top and bottom walls (zero gradient condition except at the inlet)
    inletLocIter = 1;                                                    % Used to index into the inletLocations matrix
    xPos = delta_x/2;
    for i = 2:size(pseudo_u, 2)
        xPos = delta_x/2 + (i - 2)*delta_x;
        if inletLocIter < numRecirculationInlets && xPos > inletLocations(inletLocIter, 2)
            inletLocIter = inletLocIter + 1;
        end
        if xPos > inletLocations(inletLocIter,1) && xPos < inletLocations(inletLocIter,2)
            pressure_new(end,:) = 2*p_inlet*ones(1,numel(pressure_new(end,:))) - pressure_new(end - 1,:);                       % At the top wall inlet
            pressure_new(1,:) = 2*p_inlet*ones(1,numel(pressure_new(1,:))) - pressure_new(2,:);                       % At the bottom wall inlet
        else
            % No inlet (zero gradient)
            pressure_new(end,:) = pressure_new(end - 1,:);    % Top wall
            pressure_new(1,:) = pressure_new(2,:);            % Bottom wall
        end
    end
end

