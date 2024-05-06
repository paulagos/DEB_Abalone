function captured_times = capture_EHp_times()

% Compute model predictions and plots the maturity
% at puberty (E_H >= E_Hp) for one individual at several
% temperature, food quality and pCO2 values,
% using standard DEB equations (Kooijman 2010)
%
% simu : structure with individual features and model predictions
%
% calls :   init.m
%           indiv.m
%
% 13/10/2023 - Paulo F. Lagos
%---------------------------------------------------------------
%% Load and unpack parameters
filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';
res = load(filename);
par  = res.par;          % get primary parameters
cPar = parscomp_st(par); % get compound parameters
vars_pull(par);
vars_pull(cPar);

% Define arrays for different environemntal conditions
food_array = linspace(0, 0.3, 5); % kap_X = 0.8
temp_array = linspace(12, 25, 5); % temperatues (deg C)
pco2_array = linspace(400, 1200, 5); % pCO2 (uatm)

% Initialize a data structure to store captured times
captured_times = zeros(length(food_array), length(temp_array), length(pco2_array));

% Nested loops to cover all combinations of environmental conditions
for food_idx = 1:length(food_array)
    for temp_idx = 1:length(temp_array)
        for pco2_idx = 1:length(pco2_array)
            
            % Get the current environmental conditions
            food = food_array(food_idx);
            temp = temp_array(temp_idx);
            pco2 = pco2_array(pco2_idx);
            
            % Initialize simulation with current conditions
            simu = init(filename, SSP, location, varEnv, envFoodQ, gen_number);
            
            %% 2 - calculate state variables
            simu.tEVHR(:,:,1) = indiv(simu, 1);
            
            % Check for the condition to capture and store the time
            condition_time = find_time_EHp(simu.tEVHR(:,4), par.E_Hp, simu.tEVHR(:,1));
            
            disp(condition_time);
            
            % Store the captured time for the current conditions
            captured_times(pco2_idx, temp_idx, food_idx) = condition_time;
                      
        end
    end
end

disp(captured_times);

%% code to creating plots, analyze data, etc. %%%%%%%%

% Calculate the central food quality level
central_food_quality = food_array(round(length(food_array) / 2));

% Define the food quality levels
food_levels = [min(food_array), max(food_array), central_food_quality];

% Create subplots for each food quality level
for i = 1:length(food_levels)
    food_quality_index = food_array == food_levels(i); % Choose the desired food quality level
    
    % Get a 2D slice from captured_times for the chosen food quality level
    captured_times_slice = captured_times(:, :, food_quality_index);
   
    % Create a finer grid for interpolation
    fine_temp = linspace(min(temp_array), max(temp_array), 100);
    fine_pco2 = linspace(min(pco2_array), max(pco2_array), 100);
    
    [grid_temp, grid_pco2] = meshgrid(fine_temp, fine_pco2);
    
    % Interpolate the data to the finer grid
    smoothed_captured_times = interp2(temp_array, pco2_array, captured_times_slice, grid_temp, grid_pco2, 'linear');

    % Create a 3D surface plot with the smoothed data for this food quality level
    subplot(1, 3, i);
    surf(grid_temp, grid_pco2, smoothed_captured_times);
    xlabel('Temperature (°C)');
    ylabel('pCO2 Levels (uatm)');
    zlabel('Time at puberty');
    title(['Time at puberty for Food Quality ' num2str(food_levels(i))]);
    colorbar; % Add a color bar to indicate the captured times
    
end

% Create a box plot for each slice
figure;

for i = 1:length(food_levels)
    % Extract the data for the current food quality level
    current_data = captured_times(:,:,i);
    
    % Create a box plot for the current slice
    subplot(1, 4, i); % Adjust the number of subplots as needed
    
    % Customize box plot properties
    boxplot(current_data, 'Widths', 0.5, 'OutlierSize', 5);
    
    % Set x-axis tick positions and labels to pCO2 values from pco2_array
    xtickpositions = 1:length(pco2_array);
    xticklabels(arrayfun(@num2str, pco2_array, 'UniformOutput', false));
    
    xlabel('pCO2 Levels');
    ylabel('Time at Puberty');
    title(['Food Quality Level ' num2str(food_levels(i))]);
end
    
    
% Create a box plot for each temperature level
figure;

for i = 1:length(temp_array)
    % Extract the data for the current temperature level
    current_data = squeeze(captured_times(i, :, :));
    
    % Create a box plot for the current slice
    subplot(1, 5, i); % Adjust the number of subplots as needed
    
    % Customize box plot properties
    boxplot(current_data);
    
    xlabel('Food Quality Levels');
    ylabel('Time at Puberty');
    title(['Temperature Level ' num2str(temp_array(i)) '°C']);
end
    



% Create a separate heatmap for each food level
for i = 1:length(food_levels)
    % Create a new figure for each heatmap
    figure;
    
    % Extract the data for the current food level
    captured_times_slice = captured_times(:, :, i);
    
    % Create the heatmap
    heatmap(temp_array, pco2_array, captured_times_slice);
    
    % Customize the labels and title
    xlabel('Temperature (°C)');
    ylabel('pCO2 Levels (uatm)');
    title(['Heatmap for Food Quality Level ' num2str(food_levels(i))]);
    
    % Add a color bar to indicate the captured times
    colorbar;
end



% Number of food levels
num_food_levels = length(food_levels);

% Create subplots for line plots
for i = 1:num_food_levels
    % Extract the data for the current food level
    captured_times_slice = captured_times(:, :, i);
    
    % Create a new figure for each subplot
    figure;
    
    % Create line plots for time at puberty vs. temperature
    plot(temp_array, captured_times_slice);
    
    % Customize the labels and title
    xlabel('Temperature (°C)');
    ylabel('Time at Puberty');
    title(['Time at Puberty vs. Temperature for Food Quality Level ' num2str(food_levels(i))]);
    
    % Add a legend to differentiate different pCO2 levels (if needed)
    legend(arrayfun(@num2str, pco2_array, 'UniformOutput', false), 'Location', 'Best');
end




% Create a cell array to store data for each food level
data = cell(1, 5);

% Iterate over each food level
for food_level = 1:5
    % Extract data for the current food level
    food_data = captured_times(:, :, food_level);
    
    % Convert the 2D matrix into a column vector for boxplot
    data{food_level} = food_data(:);
end

% Create a box plot
figure;
boxplot(data, 'Labels', {'Food Level 1', 'Food Level 2', 'Food Level 3', 'Food Level 4', 'Food Level 5'});
xlabel('Food Levels');
ylabel('Time at Puberty');
title('Comparison of Time at Puberty for Different Food Levels');






    


 % Calculate average and standard deviation for each temperature level
avg_time_at_puberty = mean(captured_times_slice, 3);
std_time_at_puberty = std(captured_times_slice, 0, 2);

% Create a scatter plot for the interaction between Temperature, pCO2, and Time at Puberty
figure;
hold on;

jitter = 0.3; % Adjust the jitter value as needed

for i = 1:length(temp_array)
    temp_level = temp_array(i);
    
    % Select the data for the current temperature level
    data_for_temp = captured_times_slice(i, :);
    
    % Add a jitter to the x-coordinates to reduce overlap
    jitter_values = (rand(size(data_for_temp)) - 0.5) * jitter;
    
    scatter(temp_level * ones(size(data_for_temp)) + jitter_values, data_for_temp, [], pco2_array, 'filled', 'DisplayName', ['Temperature = ' num2str(temp_level) ' °C']);
end

% Add average markers with standard deviation error bars
for i = 1:length(temp_array)
    errorbar(temp_array(i), avg_time_at_puberty(i), std_time_at_puberty(i), 'ks', 'DisplayName', ['Temperature = ' num2str(temp_array(i)) ' °C']);
end

% Customize the color of data points to represent pCO2
colormap(parula);  % You can change the colormap as needed
c = colorbar;
c.Label.String = 'pCO2 Levels (uatm)';

xlabel('Temperature (°C)');
ylabel('Time at Puberty');
title('Interaction between Temperature, pCO2, and Time at Puberty');
legend('Location', 'Northwest');


end