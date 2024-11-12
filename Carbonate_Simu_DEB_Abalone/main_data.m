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
    Food_quality = cellfun(@str2double, Food_quality);
end
if iscell(zoom_f)
    zoom_f = cellfun(@str2double, zoom_f);
end
if iscell(variable_env)
    variable_env = cellfun(@str2double, variable_env);
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
    % (other cases for other usernames as in your original script)
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

% Initialize color and legend settings
allColors = ["#77AC30", "#0072BD", "#7E2F8E", "#A2142F", "#A2143F", "#D95319", "#4DBEEE", "#EDB120", "#7E2F8E", "#FF00FF", "#00FF00", "#FF0000"];
legendLabels = {};  % Initialize as empty
colorIndex = 1;

% Define plotting function based on conditions
if length(all_locations) == 1 && length(SSP_scenarios) == 1
    plotFunc = @(simu) get_plots_simple(simu, gen_number, legendLabels, colorIndex, Food_quality, all_locations{1}, SSP_scenarios{1}, 1);
elseif length(all_locations) > 1 && length(SSP_scenarios) == 1
    plotFunc = @(simu) get_plots_medium(simu, gen_number, legendLabels, colorIndex, Food_quality, SSP_scenarios, all_locations, 1, 1, 1);
elseif length(all_locations) == 1 && length(SSP_scenarios) > 1
    plotFunc = @(simu) get_plots_medium(simu, gen_number, legendLabels, colorIndex, Food_quality, SSP_scenarios, all_locations, 1, 1, 1);
else
    plotFunc = @(simu) get_plots_complex(simu, gen_number, legendLabels, colorIndex, Food_quality, SSP_scenarios, all_locations, zoom_f, 1, 1, 1, 1);
end

% Loop over each combination of inputs for the simulation
tStart = tic;
for h = 1:length(Food_quality)
    envFoodQ = Food_quality(h);
    for i = 1:length(SSP_scenarios)
        SSP = SSP_scenarios{i};
        for k = 1:length(all_locations)
            location = all_locations{k};
            for m = 1:length(zoom_f)
                zoom = zoom_f(m);

                % Run the simulation for this combination
                simu = init(filename, SSP, location, variable_env{1}, envFoodQ, gen_number, zoom, born_date);

                % Pre-allocate various result data frames for speed
                t    = NaN(simu.tm+1, gen_number);
                L_w  = t;
                Ww_V = t;
                GSI  = t;
                
                % Continue to collect data for remaining generations
                for l = 1:gen_number
                    simu.tEVHR(:,:,m,l)   = indiv(filename, simu, l);
                    simu.tFluxes(:,:,m,l) = get_fluxes(filename, simu, m, l);
                    simu.obs(:,:,m,l)     = get_obs(simu, m, l);
                    
                    % Store data for each generation in arrays
                    t(:,l)    = simu.tEVHR(:,1,m,l);
                    L_w(:,l)  = simu.obs(:,3,m,l);
                    Ww_V(:,l) = simu.obs(:,4,m,l);
                    GSI(:,l)  = simu.obs(:,8,m,l);               
                end
                 
                t    = t(:);
                L_w  = L_w(:);
                Ww_V = Ww_V(:);
                GSI  = GSI(:);
                
                % Save results for the current scenario combination in a single file
                result_filename = sprintf('result_%s_%s_FQ%.2f_ZF%.2f_GEN%d.mat', SSP, location, envFoodQ, zoom, l);
                result_filepath = fullfile(saveDir, result_filename);
                save(result_filepath, 't', 'L_w', 'Ww_V', 'GSI');
                
                % Update legend labels and colorIndex for plotting if needed
                newLabel = sprintf('Food: %.2f', envFoodQ);
                legendLabels = unique([legendLabels, newLabel]);
                colorIndex = h;

                % Plot using selected plotting function
                plotFunc(simu);
                
            end
        end
    end
end

tEnd = toc(tStart);
fprintf('THIS SIMULATION TOOK: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));