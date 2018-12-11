function [] = plotResults(u_velocity, v_velocity, pressureField, ...
    temperature, domainLength, domainWidth,dx,dy,numRecirculationInlets)
    % Pre-processing to generate contours
%     [numRows, numCols] = size(u_velocity(:,2:end-1));
%     numNodes_x = numCols; numNodes_y = numRows;
%     x = 0:dx:domainLength;
%     y = 0:dy:domainWidth;
%     
%     % Plotting all properties at the pressure/temperature nodes
    x = dx/2:dx:domainLength;
    y = dy/2:dy:domainWidth;
    [X, Y] = meshgrid(x,y);
    
    % Calculating the velocities at the pressure nodes
    numPressureNodes_x = size(pressureField,2);             % Including ghost nodes
    numPressureNodes_y = size(pressureField,1);
    u_vel_cv = zeros(numPressureNodes_y,numPressureNodes_x);
    v_vel_cv = zeros(numPressureNodes_y,numPressureNodes_x);
    % u-velocity
%     numUVelNodes_x = size(u_velocity,2);
%     numUVelNodes_y = size(u_velocity,1);
    for j = 2:numPressureNodes_y - 1
        for i = 2:numPressureNodes_x - 1
            % Interpolating the velocities at the pressure nodes
            u_vel_cv(j,i) = 0.5*(u_velocity(j,i) + u_velocity(j,i - 1));
        end
    end
    % v-velocity
%     numVVelNodes_x = size(v_velocity,2);
%     numVVelNodes_y = size(v_velocity,1);
    for j = 2:numPressureNodes_y - 1
        for i = 2:numPressureNodes_x - 1
            % Interpolating the velocities at the pressure nodes
            v_vel_cv(j,i) = 0.5*(v_velocity(j,i) + v_velocity(j - 1,i));
        end
    end
    
    % Velocity vectors
    figure;
    u_vel_cv = u_vel_cv(2:end - 1,2:end - 1);
    v_vel_cv = v_vel_cv(2:end - 1,2:end - 1);
    q = quiver(X, Y, u_vel_cv, v_vel_cv);
%     q.AutoScale = 'on';
%     q.AutoScaleFactor = 1.5;
    hold on;
    xlabel('x', 'Interpreter', 'latex');
    ylabel('y', 'Interpreter', 'latex');
    title(strcat('Velocity vectors for \,',num2str(numRecirculationInlets),' inlets'), 'Interpreter', 'latex');
    axis([-0.5,3.5,0,1]);
    % Plotting the temperature field
%     temperature = temperature(2:end - 1,2:end - 1);
%     xScalar = dx/2:dx:domainLength;
%     yScalar = dy/2:dy:domainWidth;
%     [X_scalar,Y_scalar] = meshgrid(xScalar,yScalar);
    figure;
    [c,h] = contourf(X,Y,temperature,50);%'ShowText', 'on');
    set(h,'LineColor','none')
    xlabel('x', 'Interpreter', 'latex');
    ylabel('y', 'Interpreter', 'latex');
    title(strcat('Temperature contours for \,',num2str(numRecirculationInlets),' inlets'), 'Interpreter', 'latex');

    % Plotting the temperatures at the outlet
    figure;
    outletTemperature = temperature(:,end);
    plot(y,outletTemperature);
    title(strcat('Outlet temperature distribution for \,',num2str(numRecirculationInlets),' inlets'),'Interpreter', 'latex');
    xlabel('y','Interpreter', 'latex');
    ylabel('Outlet Temperature (C)','Interpreter', 'latex');
    ylim([250, 950]);
%     axis([0,1,0,650]);
    
%     % Streamlines
%     figure;
%     

end