% TODO - Check air physical parameters
% For numerical values - http://www.ae.iitm.ac.in/~amitk/AS1300/AS1300_TUTORIAL8_2017solution.pdf
% Solves the NS equation using primitive variable formulation and
% fractional time step method
function [] = solveCombustionChamberNS()
    maxTime = ;
    dt = ;
    % Domain parameters
    domainLength = ;
    domainWidth = ;
    leftInletVelocity = 25;
    inletTemperature = 270;                                                 % Temperature at the left inlet
    recirculation_uVel = ;                                                  % u-velocity at the recirculation inlets
    recirculation_vVel = ;                                                  % v-velocity at the recirculation inlets
    inletLocations = [];                                                    % Locations of the recirculation inlets
                                                                            % Each row represents the start and end of the inlet
    numRecirculationInlets = size(inletLocations,1);
    % Physical properties of air
    density = ;
    viscosity = ;
    specificHeat = ;
    convectionCoeff = ;
    % Discretization
    % Control volume centers correspond to pressure/temperature nodes
    % Velocity nodes are staggered wrt these nodes
    % u-velocity nodes are staggered dx/2 to the right, v-velocity nodes
    % are staggered dy/2 to the top
    numControlVols_x = ;                    
    numControlVols_y = ;
    dx = domainLength/numControlVols_x;
    dy = domainWidth/numControlVols_y;
    % Allocating the required problem variables
    pressureField = zeros(numControlVols_x + 2, numControlVols_y + 2);      % Pressure matrix includes ghost nodes too
    temperatureField = zeros(numControlVols_y, numControlVols_x);
    u_velocity = zeros(numControlVols_y + 2, numControlVols_x + 1);         % u-velocity nodes located at the left and right walls of a CV
    v_velocity = zeros(numControlVols_y + 1, numControlVols_x + 2);         % v-velocity nodes located at the bottom and top walls of a CV
    
    % Handling the boundaries
    % Setting up the inlet velocities
    u_velocity(:,1) = leftInletVelocity;
    numVNodes_x = size(v_velocity,2);
    x = dx/2;
    inletLocIter = 1;                                                       % Used to index into the inletLocations matrix
    % v-velocities at the recirculation inlets
    for i = 2:numVNodes_x - 1
        x = dx/2 + (i - 2)*dx;
        if inletLocIter < numRecirculationInlets && x > inletLocations(inletLocIter, 2)
            inletLocIter = inletLocIter + 1;
        end
        if x > inletLocations(inletLocIter,1) && x < inletLocations(inletLocIter,2)
            v_velocity(end,i) = -recirculation_vVel;                        % Top wall inlet
            v_velocity(1,i) = recirculation_vVel;                           % Top wall inlet
        end
    end
    % u-velocities at the left inlets
    u_velocity(:,1) = leftInletVelocity;
    % u-velocities at the top and bottom walls
    uVelTopWall = zeros(numControlVols_y + 2, numControlVols_x + 1);
    uVelBottomWall = uVelTopWall;
    numUNodes_x = size(u_velocity,2);
    x = 0;
    inletLocIter = 1;                                                       % Used to index into the inletLocations matrix
    for i = 1:numUNodes_x
        x = (i - 1)*dx;
        if inletLocIter < numRecirculationInlets && x > inletLocations(inletLocIter, 2)
            inletLocIter = inletLocIter + 1;
        end
        if x > inletLocations(inletLocIter,1) && x < inletLocations(inletLocIter,2)
            uVelTopWall(i) = recirculation_uVel;                            % At the top wall inlet
            uVelBottomWall(i) = recirculation_uVel;                         % At the bottom wall inlet
        end
    end
    % Temperature at the inlet
    temperatureField(:,1) = inletTemperature;
    
    currTime = 0;
    while currTime < maxTime
        % Solving for intermediate velocity field
        [intermediate_u_vel, intermediate_v_vel] = pseudo_vel_calc(u_velocity, v_velocity, ...
                                uVelBottomWall, uVelTopWall, dx, dy, dt, viscosity, density, g_x, g_y);
        % Solving for the pressure at the new time step using the
        % intermediate velocity field
        
        
        % Solving for the new velocity field using the pressure field
        [u_velocity, v_velocity] = update_vel(intermediate_u_vel, intermediate_v_vel, pressureField, dt, ...
                                                    uVelBottomWall, uVelTopWall);
                                                
        % Solving the energy equation for the temperature field
        temperatureField = solveEnergyEquation(temperatureField, u_velocity, v_velocity, pressureField,viscosity,density,specificHeat,...
                        convectionCoeff,heatSourceLocation,inletLocations,inletTemperature,dx,dy,dt);
                            
        currTime = currTime + dt;
    end
end