% TODO - Check air physical parameters
% For numerical values - http://www.ae.iitm.ac.in/~amitk/AS1300/AS1300_TUTORIAL8_2017solution.pdf
% Solves the NS equation using primitive variable formulation and
% fractional time step method
% function [] = solveCombustionChamberNS()
    clear; clc;
    maxTime = 2;
    dt = 1e-4;
    g_x = 0;
    g_y = -9.81;                                                            % gravitational constant
    % Domain parameters
    domainLength = 3;
    domainWidth = 1;
    leftInletVelocity = 0.25;
    inletTemperature = 270 + 273.13;                                                 % Temperature at the left inlet
    inletPressure = 6*101325;                                               % Inlet pressure to the chamber is 1 atm
    recirculation_uVel = 0;                                                 % u-velocity at the recirculation inlets
    recirculation_vVel = 0.25;                                               % v-velocity at the recirculation inlets
    numRecirculationInlets = 2;
    inletLocations = zeros(numRecirculationInlets,2);                       % Locations of the recirculation inlets
                                                                            % Each row represents the start and end of the inlet
    inletSize = domainLength/10;
    % Recirculation inlets uniformly distributed between domainLength/3 and
    % 4*domainLength/5
    % Each row of inletLocations stores the start and end x-coordinates of
    % the inlets
    % NEW
    inletLocations(1,:) = [0.3*domainLength, 0.3*domainLength + inletSize];
    inletLocations(2,:) = [0.5*domainLength, 0.5*domainLength + inletSize];
%     inletLocations(3,:) = [0.7*domainLength, 0.7*domainLength + inletSize];

%     % ORIGINAL
%     inletLocations(1,:) = [domainLength/3, domainLength/3 + inletSize];
%     inletLocations(2,:) = [domainLength/3 + 7*domainLength/60, domainLength/3 + 7*domainLength/60 + inletSize];
%     inletLocations(3,:) = [domainLength/3 + 2*7*domainLength/60, domainLength/3 + 2*7*domainLength/60 + inletSize];
%     inletLocations(4,:) = [domainLength/3 + 3*7*domainLength/60, domainLength/3 + 3*7*domainLength/60 + inletSize];
    % Heat source location (between 0.1*domainLength and 0.2*domainLength)
    heatSourceLocation = [0.1*domainLength, 0.3*domainLength,0.25*domainWidth,0.75*domainWidth];
    flameStart_x = heatSourceLocation(2);
    % Physical properties of air (taken at 300C)
    % source - https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
    density = 0.6158;
%     density = 1.225;
    viscosity = 2.934e-5;
    specificHeat = 1044;
    convectionCoeff = 35.45*5;                                                % source - https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
    % Reynolds number
    Re = density*leftInletVelocity*domainWidth/viscosity
    
    % source - http://bouteloup.pierre.free.fr/lica/phythe/don/air/air_k_plot.pdf
    thermalCond = @(T)(1.5207e-11*T*T*T - 4.8574e-08*T*T + 1.0184e-04*T - 3.9333e-04);

    % Discretization
    % Control volume centers correspond to pressure/temperature nodes
    % Velocity nodes are staggered wrt these nodes
    % u-velocity nodes are staggered dx/2 to the right, v-velocity nodes
    % are staggered dy/2 to the top
    numControlVols_x = 50;                    
    numControlVols_y = 20;
    dx = domainLength/numControlVols_x;
    dy = domainWidth/numControlVols_y;
    % Allocating the required problem variables
    pressureField = zeros(numControlVols_y + 2, numControlVols_x + 2);      % Pressure matrix includes ghost nodes too
    pressureField(:,1) = inletPressure;                                     % Setting the ghost nodes to the left of left wall
    pressureField(:,2) = inletPressure;                                     % and the first set of internal nodes to the inlet pressure
    temperatureField = inletTemperature*ones(numControlVols_y, numControlVols_x);
    u_velocity = leftInletVelocity*ones(numControlVols_y + 2, numControlVols_x + 1);         % u-velocity nodes located at the left and right walls of a CV
    v_velocity = zeros(numControlVols_y + 1, numControlVols_x + 2);         % v-velocity nodes located at the bottom and top walls of a CV
    
    % Handling the boundaries
    % TEMP -  to debug effect of inlet
    inletPresent = true;

    % Setting up the inlet velocities
    u_velocity(:,1) = leftInletVelocity;
    u_velocity(2:end - 1,2:end - 1) = leftInletVelocity;
    numVNodes_x = size(v_velocity,2);
    x = dx/2;
    inletLocIter = 1;                                                       % Used to index into the inletLocations matrix
    % v-velocities at the recirculation inlets
    for i = 2:numVNodes_x - 1
        x = dx/2 + (i - 2)*dx;
        if inletLocIter < numRecirculationInlets && x > inletLocations(inletLocIter, 2)
            inletLocIter = inletLocIter + 1;
        end
        if inletPresent && x > inletLocations(inletLocIter,1) && x < inletLocations(inletLocIter,2)
            v_velocity(end,i) = -recirculation_vVel;                        % Top wall inlet
            v_velocity(1,i) = recirculation_vVel;                           % Top wall inlet
%             % TEMP
%             v_velocity(end,i) = 0;                        % Top wall inlet
%             v_velocity(1,i) = 0;                           % Top wall inlet
        end
    end
    % u-velocities at the left inlets
    u_velocity(:,1) = leftInletVelocity;
    
    % u-velocities at the top and bottom walls
    uVelTopWall = zeros(1,numControlVols_x + 1);
    uVelBottomWall = uVelTopWall;
    numUNodes_x = size(u_velocity,2);
%     TEMP - CHECKING WITHOUT RECIRCULATION INLETS
    x = 0;
    inletLocIter = 1;                                                       % Used to index into the inletLocations matrix
    for i = 1:numUNodes_x
        x = (i - 1)*dx;
        if inletLocIter < numRecirculationInlets && x > inletLocations(inletLocIter, 2)
            inletLocIter = inletLocIter + 1;
        end
        if inletPresent && x > inletLocations(inletLocIter,1) && x < inletLocations(inletLocIter,2)
            uVelTopWall(i) = recirculation_uVel;                            % At the top wall inlet
            uVelBottomWall(i) = recirculation_uVel;                         % At the bottom wall inlet
        end
    end

    % Temperature at the inlet
    temperatureField(:,1) = inletTemperature;
    
    % Parameters for Gauss-Seidel with SOR
    sor_factor = 1.9397;
    epsilon = 1e-2;
    
    
    % Building the coefficient matrix to solve the Energy equation using
    % Crank-Nicolson
    K_initial = thermalCond(inletTemperature);
    coefficient = thermalCond(inletTemperature)*dt/(density*specificHeat);
    coeffMatrix = buildCoeffMatrix(coefficient,temperatureField, dx, dy);
    coeffMatInv = inv(coeffMatrix);
    currTime = dt;
    convergenceTolerance = 1e-4;
    temperatureError = 100;
    oldTemperatureField = temperatureField;
    while currTime < maxTime
%     while temperatureError > convergenceTolerance
        % Solving for intermediate velocity field
        [intermediate_u_vel, intermediate_v_vel] = pseudo_vel_calc(u_velocity, v_velocity, inletLocations,recirculation_uVel,...
                                uVelBottomWall, uVelTopWall, dx, dy, dt, viscosity, density, g_x, g_y, inletPresent);
        % Solving for the pressure at the new time step using the
        % intermediate velocity field
        pressureField = pressure_calc(pressureField, intermediate_u_vel, intermediate_v_vel,...
                            inletPressure, inletLocations, inletLocations, density, ...
                            dx, dy, dt, sor_factor, epsilon, inletPresent);
        
        % Solving for the new velocity field using the pressure field
        [u_velocity, v_velocity] = update_vel(intermediate_u_vel, intermediate_v_vel, pressureField, dt,dx,dy, ...
                                                    uVelBottomWall, uVelTopWall, density);
                                                
        % Solving the energy equation for the temperature field
        temperatureField = solveEnergyEquation(temperatureField, u_velocity, v_velocity,viscosity,density,specificHeat,...
                        convectionCoeff,thermalCond,heatSourceLocation,inletLocations,inletTemperature,coeffMatrix,dx,dy,dt, inletPresent);
                            
        disp(strcat("Current time: ",num2str(currTime)));
        % Updating the heat source location (x) as the flame spreads with time
%         % Linear traversal
%         heatSourceLocation(2) = flameStart_x + (domainLength - flameStart_x)*(currTime/maxTime);
%         % Using sigmoid
%         heatSourceLocation(2) = flameStart_x + (domainLength - flameStart_x)*(1/(1 + exp(-currTime)));
        currTime = currTime + dt;
        temperatureError = max(max(abs(temperatureField - oldTemperatureField)));
        oldTemperatureField = temperatureField;
    end
    
    plotResults(u_velocity, v_velocity, pressureField, temperatureField, domainLength,domainWidth,dx,dy);
    
% end