function [] = plotResults(u_velocity, v_velocity, pressureField, temperature, domainLength, domainWidth,dx,dy)
    % Pre-processing to generate contours
%     [numRows, numCols] = size(u_velocity(:,2:end-1));
%     numNodes_x = numCols; numNodes_y = numRows;
    x = 0:dx:domainLength;
    y = 0:dy:domainWidth;
%     x = linspace(0,domainLength,numNodes_x); y = linspace(0, domainWidth, numNodes_y);
    [X, Y] = meshgrid(x,y);
    
    % Velocity vectors
    % figure;
    u_vel = u_velocity(2:end,:);
    v_vel = v_velocity(:,2:end);
    quiver(X, Y, u_vel, v_vel);
    hold on;
    xlabel('x', 'Interpreter', 'latex');
    ylabel('y', 'Interpreter', 'latex');
    title('Velocity vectors at t = 2s', 'Interpreter', 'latex');
    
    % Plotting the temperature field
    xScalar = dx/2:dx:domainLength;
    yScalar = dy/2:dy:domainWidth;
    [X_scalar,Y_scalar] = meshgrid(xScalar,yScalar);
    figure;
    contourf(X_scalar,Y_scalar,temperature);%'ShowText', 'on');
    xlabel('x', 'Interpreter', 'latex');
    ylabel('y', 'Interpreter', 'latex');
    title('Temperatures at t = 2s', 'Interpreter', 'latex');

end