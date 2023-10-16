function captured_times = capture_condition_times()

% Compute model predictions and plots the maturity
% at metamorphosis (E_H >= E_Hj) for one individual at several
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
% filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';
% res = load(filename);
% par  = res.par;          % get primary parameters
% cPar = parscomp_st(par); % get compound parameters
% vars_pull(par);
% vars_pull(cPar);


% Prompt the user to select the .mat file
[filename, path] = uigetfile('*.mat', 'Select the .mat file');

if filename == 0
    % User canceled the selection
    disp('File selection canceled.');
    return;
end

% Construct the full file path
fullFilePath = fullfile(path, filename);

% Load the data
data = load(fullFilePath);

% Extract primary and compound parameters
par = data.par;
cPar = parscomp_st(par);

% Pull the variables into the workspace
vars_pull(par, cPar);
vars_pull(cPar);


% Define arrays for different environemntal conditions
food_array = linspace(-50, 0, 10); % kap_X
temp_array = linspace(12, 35, 10); % temperatues (deg C)
pco2_array = linspace(420, 2500, 10); % pCO2 (uatm)

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
            simu = init('RCP1', 1, food, 1, pco2, temp); % Modify init to accept additional argument
            
            %% 2 - calculate state variables
            simu.tEVHR(:,:,1) = indiv(simu, 1);
            
            % Check for the condition you want to capture and store the time
            condition_time = find_time_E_Hj(simu.tEVHR(:,4), par.E_Hj, simu.tEVHR(:,1));
            
            % Store the captured time for the current conditions
            captured_times(food_idx, temp_idx, pco2_idx) = condition_time;
            
        end
    end
end

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
    smoothed_captured_times = interp2(temp_array, pco2_array, captured_times_slice, grid_temp, grid_pco2, 'spline');
    
    % Create a 3D surface plot with the smoothed data for this food quality level
    subplot(1, 3, i);
    surf(grid_temp, grid_pco2, smoothed_captured_times);
    xlabel('Temperature (°C)');
    ylabel('pCO2 Levels (uatm)');
    zlabel('Metamorphosis Time');
    title(['Metamorphosis Time for Food Quality ' num2str(food_levels(i))]);
    colorbar; % Add a color bar to indicate the captured times
        
end


% % Create all possible combinations of grouping variables
% temp_combinations = [];
% pco2_combinations = [];
% 
% for i = 1:length(temp_array)
%     for j = 1:length(pco2_array)
%         temp_combinations = [temp_combinations; temp_array(i)];
%         pco2_combinations = [pco2_combinations; pco2_array(k)];
%     end
% end
% 
% % Perform the ANOVA with the combinations and captured_times
% [p, tbl, stats] = anovan(captured_times(:), {temp_combinations, pco2_combinations}, 'model', 'interaction', 'varnames', {'Temperature', 'pCO2'});
% 
% % Display the results of the ANOVA
% disp('Results of Three-Way ANOVA:');
% disp(tbl);

% Extract the p-values for the effects of temperature, food quality, and pCO2
% p_temp = p(1);
% p_food = p(2);
% p_pco2 = p(3);

% Box plot for metamorphosis times at different food quality levels
% figure;
% subplot(1, 2, 1);
% food_labels = {'Min', 'Max', 'Central'};
% boxplot(captured_times(:), {food_combinations_idx}, 'Labels', food_labels);
% xlabel('Food Quality Level');
% ylabel('Metamorphosis Time');
% title('Box Plot of Metamorphosis Times at Different Food Quality Levels');
% text(1, max(captured_times(:)), sprintf('p(Food) = %.4f', p_food), 'HorizontalAlignment', 'center');

% % Box plot for metamorphosis times at different pCO2 levels
% subplot(1, 2, 2);
% pco2_labels = {'Low', 'High'};
% boxplot(captured_times(:), {pco2_combinations_idx}, 'Labels', pco2_labels);
% xlabel('pCO2 Level');
% ylabel('Metamorphosis Time');
% title('Box Plot of Metamorphosis Times at Different pCO2 Levels');
% text(1, max(captured_times(:)), sprintf('p(pCO2) = %.4f', p_pco2), 'HorizontalAlignment', 'center');



% Create subplots for each food quality level
figure;

for i = 1:length(food_levels)
    food_level = food_levels(i);
    
    % Find the indices corresponding to the specific food quality level
    food_indices = find(food_array == food_level);
    
    % Extract the data for this food quality level
    captured_times_slice = captured_times(:, :, food_indices);
    temp_values = temp_array(food_indices);
    pco2_values = pco2_array(food_indices);
    
    % Create scatter plots for temperature
    subplot(3, 2, 2 * (i - 1) + 1);
    scatter(temp_values, captured_times_slice(:), 'filled');
    xlabel('Temperature (°C)');
    ylabel('Metamorphosis Time');
    title(['Food Quality Level: ' num2str(food_level) ' - Temperature']);
    
    % Create scatter plots for pCO2
    subplot(3, 2, 2 * i);
    scatter(pco2_values, captured_times_slice(:), 'filled', 'r');
    xlabel('pCO2 Levels (uatm)');
    ylabel('Metamorphosis Time');
    title(['Food Quality Level: ' num2str(food_level) ' - pCO2']);

end