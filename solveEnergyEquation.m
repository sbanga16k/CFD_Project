% Solves the energy equation for the temperature field at a particular
% timestep
% Pressures and temperatures are stored at the same nodes
% u and v - velocities are stored at a staggered location
% u velocity ghost nodes are needed above the top and below the bottom
% boundary in order to evaluate du/dy
% v velocity ghost nodes are needed before the left and beyond the right
% boundary in order to evaluate dv/dx
% heatSourceLocation = [startX, endX] for the heat source, it acts along
% the entire width of the domain as convectionCoeff*(flameTemperature - airTemperature)

% TODO - Check thermal conductivity, convection coefficient
function [temperatureField] = solveEnergyEquation(oldTemperatureField, u_velocity, v_velocity, pressure,viscosity,density,specificHeat,...
                        convectionCoeff,heatSourceLocation,dx,dy,dt)
    % The number of internal nodes in the discretization
    numInternalNodes_x = size(pressure,1);
    numInternalNodes_y = size(pressure,2);
    % Initializing the new temperature field
    temperatureField = oldTemperatureField;
    % Heat source parameters
    heatSourceStart = heatSourceLocation(1);
    heatSourceEnd = heatSourceLocation(2);
    flameTemperature = ;        % TODO
    xPos = dx/2;                                    % x-coordinate of the pressure/temperature nodes (center of CV)
    % Going over all the internal nodes to evaluate the new temperatures
    for j = 1:numInternalNodes_y
        for i = 1:numInternalNodes_x
            xPos = dx/2 + (i - 1)*dx;
            % Velocity gradients
            du_dx = (u_velocity(j + 1,i + 1) - u_velocity(j + 1,i))/dx;
            dv_dy = (v_velocity(j + 1,i + 1) - v_velocity(j,i + 1))/dx;
            % u-velocities at the top and bottom wall of a CV (interpolation)
            u_top = 0.5*(u_velocity(j + 1,i + 1) + u_velocity(j + 1,i + 1));
            u_bottom = 0.5*(u_velocity(j + 1,i + 1) + u_velocity(j,i + 1));
            % v-velocities at the left and right wall of a CV (interpolation)
            v_left = 0.5*(u_velocity(j + 1,i + 2) + u_velocity(j + 1,i + 1));
            v_right = 0.5*(u_velocity(j + 1,i + 1) + u_velocity(j + 1,i));
            du_dy = (u_top - u_bottom)/dy;
            dv_dx = (v_right - v_left)/dx;
            % Temperature gradients
            d2T_dx2 = (oldTemperatureField(j,i + 1) - 2*oldTemperatureField(j,i) + oldTemperatureField(j,i - 1))/(dx*dx);
            d2T_dy2 = (oldTemperatureField(j - 1,i) - 2*oldTemperatureField(j,i) + oldTemperatureField(j - 1,i))/(dy*dy);
            
            currTempCoeff = 1;                          % Coefficient of the current node temperature (Tji)
            % Evaluating the new temperature at (i,j)
            convectionSource = 0;
            if xPos > heatSourceStart && xPos < heatSourceEnd
               convectionSource = convectionCoeff*(flameTemperature);
               currTempCoeff = currTimeCoeff + convectionCoeff;
            end
            bodySource = 0;
            discretizedRHS = thermalCond*(d2T_dx2 + d2T_dy2) + 2*viscosity*du_dx*du_dx + viscosity*(du_dy + dv_dx)*(du_dy + dv_dx)...
                + 2*viscosity*dv_dy*dv_dy + convectionSource + bodySource;
            temperatureField(j,i) = (oldTemperatureField(j,i) + dt*discretizedRHS/(density*specificHeat))/currTempCoeff;
        end
    end
end