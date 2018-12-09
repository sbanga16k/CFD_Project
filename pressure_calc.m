% Point Gauss-Siedel method for solving pressure equation using intermediate 
% values of u,v velocities by formulating it as elliptical PDE with provision for
% Successive over-relaxation (SOR)
function [pressure_new] = pressure_calc(pressure, pseudo_u, pseudo_v,...
    p_inlet, p_top_nozzles, p_bot_nozzles,density, ...
    delta_x, delta_y, delta_t, sor_factor, epsilon)
    % Shape of pseudo_u: (n+2) x (m+1)
    % Shape of pseudo_v: (n+1) x (m+2)
    % Shape of pressure: (n+2) x (m+2)
    
    beta = delta_x/delta_y;

    % Intializing new matrix to store updated values of phi for testing
    % convergence of algorithm by computing difference with phi values at
    % previous time step
    pressure_new = zeros(size(pressure));

    % Stores max error in the matrix for testing algo convergence
    max_error = 100;

    % Iterating over all the interior nodes to update the phi values till the
    % algo converges
    while max_error > epsilon
        for j = 2:size(pseudo_v, 1)
            for i = 2:size(pseudo_u, 2)
                pressure_new(j,i) = (1 - sor_factor)* pressure(j,i) + ...
                    (0.5*sor_factor/(1 + beta*beta))*(pressure_new(j, i-1) + ...
                    pressure(j, i+1) + ...
                    (beta*beta)*(pressure_new(j-1, i) + pressure(j+1,i))) - ...
                    0.5/(1 + beta*beta) * (density*delta_x^2/delta_t) * ...
                    ((pseudo_u(j,i) - pseudo_u(j,i-1))/delta_x + ...
                    (pseudo_v(j,i) - pseudo_v(j-1,i))/delta_y);
            end
        end
        
        % % Ghost nodes update for pressure
       
        % Left wall Ghost nodes update
        % to have the average to be equal to the value of pressure at inlet (p_inlet)
        pressure_new(:,1) = 2*p_inlet - pressure_new(:,2);

        % Right wall Ghost nodes update (assuming zero flux)
        pressure_new(:,end) = pressure_new(:,end-1);

        % Bottom wall Ghost nodes update (assuming zero flux)
        pressure_new(1,:) = pressure_new(2,:);
        for k = 1:size(p_bot_nozzles, 1)
            nozzle_inds = p_bot_nozzles(k,:)./delta_x + 1;
            pressure_new(1,nozzle_inds(1):nozzle_inds(2)) = 2*p_inlet - ...
                pressure_new(2,nozzle_inds(1):nozzle_inds(2));
        end

        % Top wall Ghost nodes update (assuming zero flux)
        pressure_new(end,:) = pressure_new(end-1,:);
        for l = 1:size(p_top_nozzles, 1)
            nozzle_inds = p_top_nozzles(l,:)./delta_x + 1;
            pressure_new(end,nozzle_inds(1):nozzle_inds(2)) = 2*p_inlet - ...
                pressure_new(end-1,nozzle_inds(1):nozzle_inds(2));
        end

        % Compute error in values between consecutive time steps and sets
        % pressure to be equal to pressure_new
        max_error = max(abs(pressure_new(:) - pressure(:)));
        pressure = pressure_new;
    end
    
end

