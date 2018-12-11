% Combustion source - Industrial Gas Turbines, Razak (Chapter 6)
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
% TODO - Ghost nodes for temperature
function [temperatureField] = solveEnergyEquation(oldTemperatureField, u_velocity, v_velocity,viscosity,density,specificHeat,...
                        convectionCoeff,thermalCond,heatSourceLocation,inletLocations,inletTemperature,coeffMatrix,dx,dy,dt,inletPresent)
    % The number of internal nodes in the discretization
    numInternalNodes_x = size(u_velocity,2) - 1;
    numInternalNodes_y = size(v_velocity,1) - 1;
    % Initializing the new temperature field
    temperatureField = oldTemperatureField;
    % Heat source parameters
    heatSourceStart_x = heatSourceLocation(1);
    heatSourceEnd_x = heatSourceLocation(2);
    heatSourceStart_y = heatSourceLocation(3);
    heatSourceEnd_y = heatSourceLocation(4);
    
    flameTemperature = 2200;
%    flameTemperature = 5100; 
    % Total number of recirculation inlets
    numRecirculationInlets = size(inletLocations,1);
% ORIGINAL
    % Going over all the internal nodes to evaluate the new temperatures
    for j = 1:numInternalNodes_y
        inletLocIter = 1;                                                       % Used to index into the inletLocations matrix
        yPos = dy/2 + (j - 1)*dy;                                               % y-coordinate of the pressure/temperature nodes (center of CV)
        for i = 1:numInternalNodes_x
            xPos = dx/2 + (i - 1)*dx;                                           % x-coordinate of the pressure/temperature nodes (center of CV)
            if inletLocIter < numRecirculationInlets && xPos > inletLocations(inletLocIter, 2)
                inletLocIter = inletLocIter + 1;
            end
            % Velocity gradients
            du_dx = (u_velocity(j + 1,i + 1) - u_velocity(j + 1,i))/dx;
            dv_dy = (v_velocity(j + 1,i + 1) - v_velocity(j,i + 1))/dx;
            
            
            % u-velocities at the top and bottom wall of a CV (interpolation)
%             u_top = 0.5*(u_velocity(j + 1,i + 1) + u_velocity(j + 1,i + 1));
%             u_bottom = 0.5*(u_velocity(j + 1,i + 1) + u_velocity(j,i + 1));
            u_top = 0.25*(u_velocity(j + 1,i + 1) + u_velocity(j + 1,i) + u_velocity(j + 2,i) + u_velocity(j + 2,i + 1));
            u_bottom = 0.25*(u_velocity(j + 1,i + 1) + u_velocity(j,i + 1) + u_velocity(j,i) + u_velocity(j + 1,i));
            
            % v-velocities at the left and right wall of a CV (interpolation)
%             v_left = 0.5*(v_velocity(j + 1,i + 2) + v_velocity(j + 1,i + 1));
%             v_right = 0.5*(v_velocity(j + 1,i + 1) + v_velocity(j + 1,i));
            v_left = 0.25*(v_velocity(j + 1,i + 1) + v_velocity(j,i + 1) + v_velocity(j,i) + v_velocity(j,i + 1));
            v_right = 0.25*(v_velocity(j + 1,i + 1) + v_velocity(j + 1,i + 2) + v_velocity(j,i + 2) + v_velocity(j,i + 1));
            du_dy = (u_top - u_bottom)/dy;
            dv_dx = (v_right - v_left)/dx;
            % Temperature gradients
            % Bottom wall
            if j == 1
                T_bottom = oldTemperatureField(j,i);
                % Checking for an inlet
                if inletPresent && xPos > inletLocations(inletLocIter,1) && xPos < inletLocations(inletLocIter,2)
                    T_bottom = 2*inletTemperature - oldTemperatureField(j,i);
                end
            else
                T_bottom = oldTemperatureField(j - 1,i);
            end
            % Top wall
            if j == numInternalNodes_y
                T_top = oldTemperatureField(j,i);
                % Checking for an inlet
                if inletPresent && xPos > inletLocations(inletLocIter,1) && xPos < inletLocations(inletLocIter,2)
                    T_top = 2*inletTemperature - oldTemperatureField(j,i);
                end
            else
                T_top = oldTemperatureField(j + 1,i);
            end
            % Left wall (inlet)
            if i == 1
                T_left = 2*inletTemperature - oldTemperatureField(j,i);
            else
                T_left = oldTemperatureField(j,i - 1);
            end
            % Right wall (zero gradient condition)
            if i == numInternalNodes_x
                T_right = oldTemperatureField(j,i);
            else
                T_right = oldTemperatureField(j,i + 1);
            end
            d2T_dx2 = (T_right - 2*oldTemperatureField(j,i) + T_left)/(dx*dx);
            d2T_dy2 = (T_top - 2*oldTemperatureField(j,i) + T_bottom)/(dy*dy);
            
            currTempCoeff = 1;                          % Coefficient of the current node temperature (Tji)
            % Evaluating the new temperature at (i,j)
            convectionSource = 0;
            % ORIGINAL Heat source formunlation
%             if xPos > heatSourceStart_x && xPos < heatSourceEnd_x && ...
%                     yPos > heatSourceStart_y && yPos < heatSourceEnd_y
%                convectionSource = convectionCoeff*(flameTemperature);
%                currTempCoeff = currTempCoeff + convectionCoeff*dt/(density*specificHeat);
%             end

            % NEW heat source formulation
            % Source dies down along the length
%             heatSourceStart_y = heatSourceLocation(3)
%             heatSourceEnd_y = heatSourceLocation(4);
            if xPos > heatSourceStart_x && yPos > heatSourceStart_y && yPos < heatSourceEnd_y
                % Using a linear decay of the heat source
%                 convectionSource = convectionCoeff*(heatSourceStart_x*flameTemperature/(xPos));
                % Using an exponential decay for the heat source
                convectionSource = convectionCoeff*(flameTemperature/(1 + exp(-heatSourceStart_x/xPos)));
%                 convectionSource = convectionCoeff*(flameTemperature);
%                 currTempCoeff = currTempCoeff + convectionCoeff*dt/(density*specificHeat);
            end
            bodySource = -v_velocity(j + 1,i + 1)*(9.81);
            K = 1000*thermalCond(oldTemperatureField(j,i));
            discretizedRHS = K*(d2T_dx2 + d2T_dy2) + 2*viscosity*du_dx*du_dx + viscosity*(du_dy + dv_dx)*(du_dy + dv_dx)...
                + 2*viscosity*dv_dy*dv_dy + convectionSource + bodySource;
%             % TEMP
%             discretizedRHS = K*(d2T_dx2 + d2T_dy2) + ...
%                 convectionSource + bodySource;

            temperatureField(j,i) = (oldTemperatureField(j,i) + dt*discretizedRHS/(density*specificHeat))/currTempCoeff;
        end
    end
        

end