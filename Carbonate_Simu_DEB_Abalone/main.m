%---------------------------------------------------------------
% Compute model predictions for the European abalone (Haliotis tuberculata)
% using parameter estimates for the European abalone.
%---------------------------------------------------------------

% Read simulation options from the MATLAB base workspace
SSP_scenarios = evalin('base', 'selectedScenarios');
all_locations = evalin('base', 'selectedLocations');
Food_quality = evalin('base', 'selectedFoodQuality');
variable_env = evalin('base', 'selectedEnvVariability');
zoom_f = evalin('base', 'zoomFactor');
gen_number = evalin('base', 'genNumber');
born_date = evalin('base', 'bornDate');
% Convert from cell of strings to numeric array (if needed)
if iscell(Food_quality)
    Food_quality = cellfun(@str2double, Food_quality);  % Convert all values to numbers
end

if iscell(zoom_f)
    zoom_f = cellfun(@str2double, zoom_f);  % Convert all values to numbers
end

if iscell(variable_env)
    variable_env = cellfun(@str2double, variable_env);  % Convert all values to numbers
end


% Ensure single selections are handled like multiple options (wrap in cell arrays if necessary)
if ~iscell(SSP_scenarios), SSP_scenarios = {SSP_scenarios}; end
if ~iscell(all_locations), all_locations = {all_locations}; end
if ~iscell(variable_env), variable_env = {variable_env}; end

% Set paths for saving results (adjust paths accordingly)
username = 'paulo';
switch username
    case 'paulo'
        filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';
        results  = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/Simulation_DEB_Abalone/Simulation_files/DEB_model v.1.0/Results/';
    case 'laure'
        filename = 'to be defined by Laure';
        results  = 'to be defined by Laure';
    case 'jonathan'
        filename = 'to be defined by Jonathan';
        results  = 'to be defined by Jonathan';
    case 'sabine'
        filename = 'to be defined by Sabine';
        results  = 'to be defined by Sabine';
    otherwise
        error('Unknown username. Please define filename and results paths.');
end

% Check if filename exists
if ~exist(filename, 'file')
    error('The environmental data file specified by filename does not exist. Please check the path.');
end

timeStamp = datestr(now, 'yyyymmmdd');
saveDir = [results, 'simu_abal_', timeStamp, '/'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

doPlotsSet1 = 1; % Set to 1 to plot SET 1, set to 0 to disable.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Set simulation Scenarios  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
% Number of total scenarios (for progress tracking)
numberOfScenarios = length(SSP_scenarios) * length(all_locations) * ...
    length(Food_quality) * length(variable_env) * length(zoom_f);
scenNo = 1;
close all;

% Initialize the legend labels (to be updated dynamically for each simulation)
legendLabels = {};  % Initialize as empty

% Initialize colorIndex to assign a unique color for each simulation
colorIndex = 1;

% Predefine a larger set of colors if necessary to avoid repetition
allColors = ["#77AC30", "#0072BD", "#7E2F8E", "#A2142F", "#A2143F", "#D95319", "#4DBEEE", "#EDB120", "#7E2F8E", "#FF00FF", "#00FF00", "#FF0000"];

% Initialize storage for simulation results
all_lengths = []; % To store all lengths for different locations and SSPs
all_weights = [];
all_GSI = [];
group_labels = {}; % To store labels for each group (e.g., Location + SSP combination)

% Loop over each combination of inputs for the simulation

for h = 1:length(Food_quality)
    envFoodQ = Food_quality(h);  % Now, this will be numeric
    
    % Run the simulation for each combination of SSP and location
    for i = 1:length(SSP_scenarios)
        SSP = SSP_scenarios{i};
        
        for k = 1:length(all_locations)
            location = all_locations{k};
            
            for m = 1:length(zoom_f)
                zoom = zoom_f(m);
                
                % Run simulation for this combination
                simu = init(filename, SSP, location, variable_env{1}, envFoodQ, gen_number, zoom, born_date); % filename passed here
      
                for l = 1:gen_number
                    simu.tEVHR(:,:,m,l) = indiv(filename, simu, l);
                    simu.tFluxes(:,:,m,l) = get_fluxes(filename, simu, m, l);
                    simu.obs(:,:,m,l) = get_obs(simu, m, l);
              
                end
                
                % Collect results for statistical analysis
                t = simu.tEVHR(:, 1);  % Time vector
                L_w = simu.obs(:, 3);  % Physical length
                Ww_V = simu.obs(:, 4);  % Wet-weight
                GSI = simu.obs(:, 8);  % GSI
                
                % Append simulation results for statistical analysis
                all_lengths = [all_lengths; L_w];
                all_weights = [all_weights; Ww_V];
                all_GSI = [all_GSI; GSI];
                
                % Add a label for this group, for statistical analysis
                group_labels = [group_labels; repmat({sprintf('Loc: %s, SSP: %s', location, SSP)}, length(L_w), 1)];
                
                % Determine which plotting function to use based on user selection
                
                % Case 1: Single SSP and single location
                if length(all_locations) == 1 && length(SSP_scenarios) == 1
                    
                    % Add legend label for this specific combination
                    newLabel     = sprintf('Food: %.2f', envFoodQ);
                    legendLabels = unique([legendLabels,newLabel]);  % Append to legendLabels
                    colorIndex   = h;  % Update color index for the next plot
                    
                    get_plots_simple(simu, gen_number, legendLabels, colorIndex, Food_quality, location, SSP, h);
                    
                    % Case 2: Multiple locations, one SSP scenario
                elseif length(all_locations) > 1 && length(SSP_scenarios) == 1
                    
                    % Add legend label for this specific combination
                    newLabel     = sprintf('Food: %.2f', envFoodQ);
                    legendLabels = unique([legendLabels,newLabel]);  % Append to legendLabels
                    colorIndex   = h;  % Update color index for the next plot
                    
                    
                    get_plots_medium(simu, gen_number, legendLabels, colorIndex, Food_quality, SSP_scenarios, all_locations, h, i, k);
                    
                    % Case 3: One location, multiple SSP scenarios
                elseif length(all_locations) == 1 && length(SSP_scenarios) > 1
                    
                    % Add legend label for this specific combination
                    newLabel     = sprintf('Food: %.2f', envFoodQ);
                    legendLabels = unique([legendLabels,newLabel]);  % Append to legendLabels
                    colorIndex   = h;  % Update color index for the next plot
                    
                    
                    get_plots_medium(simu, gen_number, legendLabels, colorIndex, Food_quality, SSP_scenarios, all_locations, h, i, k);
                    
                    % Case 4: Multiple locations and multiple SSP scenarios
                else
                    % Add legend label for this specific combination
                    newLabel     = sprintf('Food: %.2f', envFoodQ);
                    legendLabels = unique([legendLabels,newLabel]);  % Append to legendLabels
                    colorIndex   = h;  % Update color index for the next plot
                    
                    get_plots_complex(simu, gen_number, legendLabels, colorIndex, Food_quality, SSP_scenarios, all_locations, zoom_f, h, i, k, m);
                end
                
            end
        end
    end
end


% Save simulation results for statistical analysis
save('simulation_results.mat', 'all_lengths', 'all_weights', 'all_GSI', 'group_labels');

toc;

