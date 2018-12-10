%  Building a sparse matrix for solving linear system, Ax = b (Energy equation)
%  Generates the coefficient matrix(A)
%  Flags are used to check if flux is specified at any boundary
%  coefficient = dt/(density*specificHeat)
% TODO - CHECK MATRIX INVERTIBILITY
function [coeffMatrix] = buildCoeffMatrix(coefficient,temperature, dx, dy)
    numNodes_x = size(temperature,2);
    numNodes_y = size(temperature,1);
    totalNodes = numNodes_x*numNodes_y;
    coeffMatrix = zeros(totalNodes,totalNodes);

% 	 Elements are filled in the matrix in a row-major order
% 	 starting from y = 0 to y = Y_MAX
    x_coeff = -coefficient/(2*dx*dx); y_coeff = -coefficient / (2*dy*dy);
    curr_coeff = 1.0 - x_coeff - y_coeff;
    rowIdx = 1;
    for yIdx = 1:numNodes_y
        for xIdx = 1:numNodes_x
            currIdx = xIdx + (yIdx - 1)*numNodes_x;
            % T(i, j), LEFT
            coeffMatrix(rowIdx,currIdx) = coeffMatrix(rowIdx,currIdx) + curr_coeff;
            % T(i - 1, j)
            if xIdx == 1
                % Left boundary is inlet
                coeffMatrix(rowIdx,currIdx) = coeffMatrix(rowIdx,currIdx) - x_coeff;
            else
                coeffMatrix(rowIdx,currIdx) = coeffMatrix(rowIdx,currIdx - 1) + x_coeff;
            end
            % T(i + 1, j), RIGHT
            if xIdx == numNodes_x
                % Zero gradient condition at right outlet
                coeffMatrix(rowIdx,currIdx) = coeffMatrix(rowIdx,currIdx) + x_coeff;
            else
                coeffMatrix(rowIdx,currIdx) = coeffMatrix(rowIdx,currIdx + 1) + curr_coeff;
            end
            % u(i, j - 1), BOTTOM (TODO - inlet)
            if yIdx == 1
                % Zero gradient condition at right outlet
                coeffMatrix(rowIdx,currIdx) = coeffMatrix(rowIdx,currIdx) + y_coeff;
            else
                coeffMatrix(rowIdx,currIdx) = coeffMatrix(rowIdx,currIdx - numNodes_x) + y_coeff;
            end
            % u(i, j + 1), TOP (TODO - inlet)
            if yIdx == numNodes_y
                % Zero gradient condition at right outlet
                coeffMatrix(rowIdx,currIdx) = coeffMatrix(rowIdx,currIdx) + y_coeff;
            else
                coeffMatrix(rowIdx,currIdx) = coeffMatrix(rowIdx,currIdx + numNodes_x) + y_coeff;
            end
        end
        rowIdx = rowIdx + 1;
    end
end

