% Updates the state of every non-boundary cell.
% All boundary cells are of 
%function [] = cell_sim(square_dim, exp_length)

    % Create matrix of dim x dim
    square_dim = 35;
    exp_length = 20;
    center_coord = square_dim / 2;
    
    patterned_cell = cell(square_dim);
    
    % Temporary initial cell values.
    boundary_cell = [0 0 0 0 0 0 0 0 0 0];
    radius = sqrt(square_dim/2);
    
    % Build initial matrix
    for i=1:square_dim
        for j=1:square_dim
            if (sqrt((i-center_coord)^2 + (j-center_coord)^2)<(radius-.5)^2)
                patterned_cell{i,j}=initial_cell(i, j, square_dim);
            else
                patterned_cell{i,j}=boundary_cell;
            end
        end
    end
    
    % Visualize Initial State
    
    % After creation of inital matrix begin time_ticks
    
    data_arr = cell(1,exp_length);
    running_arr = patterned_cell;
    
    for i=1:exp_length
        % Tick all cells
        running_arr = tick_cell_indicators(running_arr, square_dim);
        data_arr{i} = running_arr;
   end
    
%end