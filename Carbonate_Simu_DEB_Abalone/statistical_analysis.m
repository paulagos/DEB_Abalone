
% Directory where simulation result files are saved
resultDir = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/Simulation_DEB_Abalone/Simulation_files/DEB_model v.1.0/Results/simu_abal_2024Nov06'; 
resultFiles = dir(fullfile(resultDir, 'result_*.mat'));

% Initialize arrays to store aggregated data
all_lengths = [];
all_weights = [];
all_GSI = [];
group_labels = {};

% Loop through each result file and load its data
for file = resultFiles'
    % Load the data from each result file
    data = load(fullfile(resultDir, file.name), 't', 'L_w', 'Ww_V', 'GSI');
    % Append data to the main arrays
    all_lengths = [all_lengths; data.L_w];
    all_weights = [all_weights; data.Ww_V];
    all_GSI = [all_GSI; data.GSI];
    
    % Extract group information from filename
    [~, filename, ~] = fileparts(file.name);
    tokens = regexp(filename, 'result_(?<SSP>.*?)_(?<Location>.*?)_FQ(?<FQ>[\d.]+)_ZF(?<ZF>[\d.]+)_GEN(?<GEN>\d+)', 'names');
    if ~isempty(tokens)
        group_label = sprintf('Loc: %s, SSP: %s, FQ: %s, ZF: %s, Gen: %s', ...
            tokens.Location, tokens.SSP, tokens.FQ, tokens.ZF, tokens.GEN);
        group_labels = [group_labels; repmat({group_label}, size(data.L_w, 1), 1)];
    end
     
end

% Save aggregated data for analysis
save('aggregated_simulation_results.mat', 'all_lengths', 'all_weights', 'all_GSI', 'group_labels');

% Create a table from the data
aggregated_data_table = table(group_labels, all_lengths, all_weights, all_GSI, ...
    'VariableNames', {'Group', 'Length', 'Weight', 'GSI'});

% Save as CSV
writetable(aggregated_data_table, 'aggregated_simulation_results.csv');

% Load the simulation results
load('aggregated_simulation_results.mat', 'all_lengths', 'all_weights', 'all_GSI', 'group_labels');

% Display group_labels for inspection
disp('Contents of group_labels:');
disp(group_labels);  % Helps you inspect the data format

% Initialize cell arrays for locations and scenarios
locations = cell(size(group_labels));
scenarios = cell(size(group_labels));

% Loop through group_labels to extract locations and scenarios
for i = 1:length(group_labels)
    % Extract location (text after 'Loc: ' and before ',')
    loc_start = strfind(group_labels{i}, 'Loc: ') + length('Loc: ');
    loc_end = strfind(group_labels{i}, ',') - 1;
    locations{i} = strtrim(group_labels{i}(loc_start:loc_end));
    
    % Extract scenario (text after 'SSP: ')
    ssp_start = strfind(group_labels{i}, 'SSP: ') + length('SSP: ');
    scenarios{i} = strtrim(group_labels{i}(ssp_start:end));
end

% Convert to categorical for better handling in ANOVA
locations = categorical(locations);
scenarios = categorical(scenarios);

% Combine location and scenario for group labels
group_labels_combined = strcat(string(locations), ' - ', string(scenarios));

% ------------------- Two-Way ANOVA for Length with anovan --------------------
fprintf('Performing Two-Way ANOVA for Lengths with anovan...\n');
[p_length, tbl_length, stats_length] = anovan(all_lengths, {locations, scenarios}, 'model', 'interaction', 'varnames', {'Location', 'Scenario'});
disp('ANOVA table for Length:');
disp(tbl_length);

% Plot ANOVA Table for Length
figure;
uitable('Data', tbl_length, 'ColumnName', tbl_length(1,:), 'RowName', [], 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
title('ANOVA Table for Length', 'FontSize', 14, 'FontWeight', 'bold');

% Tukey HSD post-hoc test for Length
fprintf('Performing Tukey HSD for Length...\n');
[c_length,~,~,gnames_length] = multcompare(stats_length, 'CType', 'tukey-kramer');
% Plot results of Tukey HSD for Length
figure;
multcompare(stats_length, 'CType', 'tukey-kramer');
title('Tukey HSD for Length', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;
% Save the Tukey HSD figure for Length
saveas(gcf, 'TukeyHSD_Length_NatureComm.tiff', 'tiff');

% ------------------- Two-Way ANOVA for Weight with anovan --------------------
fprintf('Performing Two-Way ANOVA for Weight with anovan...\n');
[p_weight, tbl_weight, stats_weight] = anovan(all_weights, {locations, scenarios}, 'model', 'interaction', 'varnames', {'Location', 'Scenario'});
disp('ANOVA table for Weight:');
disp(tbl_weight);

% Plot ANOVA Table for Weight
figure;
uitable('Data', tbl_weight, 'ColumnName', tbl_weight(1,:), 'RowName', [], 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
title('ANOVA Table for Weight', 'FontSize', 14, 'FontWeight', 'bold');

% Tukey HSD post-hoc test for Weight
fprintf('Performing Tukey HSD for Weight...\n');
[c_weight,~,~,gnames_weight] = multcompare(stats_weight, 'CType', 'tukey-kramer');
% Plot results of Tukey HSD for Weight
figure;
multcompare(stats_weight, 'CType', 'tukey-kramer');
title('Tukey HSD for Weight', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;
% Save the Tukey HSD figure for Weight
saveas(gcf, 'TukeyHSD_Weight_NatureComm.tiff', 'tiff');

% ------------------- Two-Way ANOVA for GSI with anovan --------------------
fprintf('Performing Two-Way ANOVA for GSI with anovan...\n');
[p_GSI, tbl_GSI, stats_GSI] = anovan(all_GSI, {locations, scenarios}, 'model', 'interaction', 'varnames', {'Location', 'Scenario'});
disp('ANOVA table for GSI:');
disp(tbl_GSI);

% Plot ANOVA Table for GSI
figure;
uitable('Data', tbl_GSI, 'ColumnName', tbl_GSI(1,:), 'RowName', [], 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
title('ANOVA Table for GSI', 'FontSize', 14, 'FontWeight', 'bold');

% Tukey HSD post-hoc test for GSI
fprintf('Performing Tukey HSD for GSI...\n');
[c_GSI,~,~,gnames_GSI] = multcompare(stats_GSI, 'CType', 'tukey-kramer');
% Plot results of Tukey HSD for GSI
figure;
multcompare(stats_GSI, 'CType', 'tukey-kramer');
title('Tukey HSD for GSI', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;
% Save the Tukey HSD figure for GSI
saveas(gcf, 'TukeyHSD_GSI_NatureComm.tiff', 'tiff');

% ------------------- Jittered Scatter Plot with Translucency for Length --------------------
figure;
scatter(categorical(group_labels_combined), all_lengths, 100, 'filled', 'jitter', 'on', 'jitterAmount', 0.2, 'MarkerFaceAlpha', 0.5);
title('Jittered Scatter Plot of Length across Locations and SSP Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Length (cm)', 'FontSize', 12);
xlabel('Group (Location + SSP)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;


% ------------------- Jittered Scatter Plot with Translucency for Weight --------------------
figure;
scatter(categorical(group_labels_combined), all_weights, 100, 'filled', 'jitter', 'on', 'jitterAmount', 0.2, 'MarkerFaceAlpha', 0.5);
title('Jittered Scatter Plot of Weight across Locations and SSP Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Weight (g)', 'FontSize', 12);
xlabel('Group (Location + SSP)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;


% ------------------- Jittered Scatter Plot with Translucency for GSI --------------------
figure;
scatter(categorical(group_labels_combined), all_GSI, 100, 'filled', 'jitter', 'on', 'jitterAmount', 0.2, 'MarkerFaceAlpha', 0.5);
title('Jittered Scatter Plot of GSI across Locations and SSP Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('GSI', 'FontSize', 12);
xlabel('Group (Location + SSP)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;



% ------------------- Customized Violin Plots --------------------

% ------------------- Customized Violin Plots for Length, Weight, and GSI --------------------

% Assuming you already have the 'violinplot' function from Daniel Ellis

% Combine location and scenario for group labels
group_labels_combined = strcat(string(locations), ' - ', string(scenarios));

% ------------------- Violin Plot for Length --------------------
figure;
vp = violinplot(all_lengths, group_labels_combined);  % Violin plot for Length
title('Customized Violin Plot of Length across Locations and SSP Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Length (cm)', 'FontSize', 12);
xlabel('Group (Location + SSP)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;

% Customizing each violin in the plot
for i = 1:length(vp)
    % Reduce size and add alpha transparency to scatter points
    vp(i).ScatterPlot.SizeData = 20;   % Smaller dots
    vp(i).ScatterPlot.MarkerFaceAlpha = 0.5;  % 50% transparency
    
    % Customize the median marker (make larger, full red)
    vp(i).MedianPlot.SizeData = 100;   % Larger mean marker
    vp(i).MedianPlot.MarkerFaceColor = 'r';  % Red color for the mean
    vp(i).MedianPlot.MarkerEdgeColor = 'r';  % Red color for mean dot edge
    
    % Customize deviation lines (thinner, black)
    vp(i).WhiskerPlot.LineWidth = 1;   % Thinner whisker (deviation) lines
    vp(i).WhiskerPlot.Color = 'k';     % Black whisker lines
    
    % Thinner and black IQR line (wider gray line)
    vp(i).BoxPlot.LineWidth = 1;       % Thinner IQR line
    vp(i).BoxPlot.Color = 'k';         % Black color for IQR box line
end


% ------------------- Violin Plot for Weight --------------------
figure;
vp = violinplot(all_weights, group_labels_combined);  % Violin plot for Weight
title('Customized Violin Plot of Weight across Locations and SSP Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Weight (g)', 'FontSize', 12);
xlabel('Group (Location + SSP)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;

for i = 1:length(vp)
    % Reduce size and add alpha transparency to scatter points
    vp(i).ScatterPlot.SizeData = 20;   % Smaller dots
    vp(i).ScatterPlot.MarkerFaceAlpha = 0.5;  % 50% transparency
    
    % Customize the median marker (make larger, full red)
    vp(i).MedianPlot.SizeData = 100;   % Larger mean marker
    vp(i).MedianPlot.MarkerFaceColor = 'r';  % Red color for the mean
    vp(i).MedianPlot.MarkerEdgeColor = 'r';  % Red color for mean dot edge
    
    % Customize deviation lines (thinner, black)
    vp(i).WhiskerPlot.LineWidth = 1;   % Thinner whisker (deviation) lines
    vp(i).WhiskerPlot.Color = 'k';     % Black whisker lines
    
    % Thinner and black IQR line (wider gray line)
    vp(i).BoxPlot.LineWidth = 1;       % Thinner IQR line
    vp(i).BoxPlot.Color = 'k';         % Black color for IQR box line
end



% ------------------- Violin Plot for GSI --------------------
figure;
vp = violinplot(all_GSI, group_labels_combined, ...
    'MarkerSize', 8, ...          % Smaller scatter points
    'ViolinAlpha', 0.1, ...        % 50% transparency
    'ViolinColor', locationColors, ... % Different colors for each location
    'MedianMarkerSize', 40, ...   % Larger median marker
    'ViolinColor', [0.8, 0.8, 0.8], ... % Light gray violin fill
    'EdgeColor', [0, 0, 0], ...          % Black edge for violins
    'BoxColor', [0, 0, 0], ...           % Black whisker and box lines
    'ShowData', true, ...          % Show data points
    'ShowMedian', true,...           % Show median marker
    'ShowWhiskers', true,...
    'ShowNotches', false);
title('Customized Violin Plot of GSI across Locations and SSP Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('GSI', 'FontSize', 12);
xlabel('Group (Location + SSP)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;


% ------------------- Save Figures --------------------

% Save the boxplots and violin plots in high resolution suitable for publication
saveas(gcf, 'Boxplot_GSI_NatureComm.tiff', 'tiff');
saveas(gcf, 'Boxplot_Length_NatureComm.tiff', 'tiff');
saveas(gcf, 'Boxplot_Weight_NatureComm.tiff', 'tiff');
saveas(gcf, 'JitteredPlot_GSI_Translucent_NatureComm.tiff', 'tiff');
saveas(gcf, 'JitteredPlot_Weight_Translucent_NatureComm.tiff', 'tiff');
saveas(gcf, 'JitteredPlot_Length_Translucent_NatureComm.tiff', 'tiff');
saveas(gcf, 'ViolinPlot_Weight_NatureComm.tiff', 'tiff');
saveas(gcf, 'ViolinPlot_Length_NatureComm.tiff', 'tiff');
saveas(gcf, 'Customized_ViolinPlot_Length_NatureComm.tiff', 'tiff');
saveas(gcf, 'Customized_ViolinPlot_Weight_NatureComm.tiff', 'tiff');
saveas(gcf, 'Customized_ViolinPlot_GSI_NatureComm.tiff', 'tiff');




% Define the number of unique combinations of locations and SSP scenarios
uniqueGroups = unique(group_labels_combined);
numGroups = length(uniqueGroups);

% Generate distinct colors for each unique combination
groupColors = num2cell(lines(numGroups), 2);  % Using 'lines' colormap and splitting into cell array

% ------------------- Violin Plot for Length --------------------
figure;
vp = violinplot(all_lengths, group_labels_combined, ...
    'MarkerSize', 20, ...          % Smaller scatter points
    'ViolinAlpha', 0.5, ...        % 50% transparency
    'MedianMarkerSize', 100, ...   % Larger median marker
    'ViolinColor', groupColors, ... % Different colors for each unique group (location + SSP)
    'EdgeColor', [0, 0, 0], ...    % Black edge for violins (RGB triplet)
    'BoxColor', [0, 0, 0], ...     % Black whisker and box lines (RGB triplet)
    'ShowData', true, ...          % Show data points
    'ShowMedian', true);           % Show median marker
title('Customized Violin Plot of Length across Locations and SSP Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Length (cm)', 'FontSize', 12);
xlabel('Group (Location + SSP)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;



% ------------------- Violin Plot for Weight --------------------
figure;
vp = violinplot(all_weights, group_labels_combined, ...
    'MarkerSize', 20, ...          % Smaller scatter points
    'ViolinAlpha', 0.5, ...        % 50% transparency
    'MedianMarkerSize', 100, ...   % Larger median marker
    'ViolinColor', groupColors, ... % Different colors for each unique group (location + SSP)
    'EdgeColor', [0, 0, 0], ...    % Black edge for violins (RGB triplet)
    'BoxColor', [0, 0, 0], ...     % Black whisker and box lines (RGB triplet)
    'ShowData', true, ...          % Show data points
    'ShowMedian', true);           % Show median marker
title('Customized Violin Plot of Weight across Locations and SSP Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Weight (g)', 'FontSize', 12);
xlabel('Group (Location + SSP)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;



% ------------------- Violin Plot for GSI --------------------
figure;
vp = violinplot(all_GSI, group_labels_combined, ...
    'MarkerSize', 20, ...          % Smaller scatter points
    'ViolinAlpha', 0.5, ...        % 50% transparency
    'MedianMarkerSize', 100, ...   % Larger median marker
    'ViolinColor', groupColors, ... % Different colors for each unique group (location + SSP)
    'EdgeColor', [0, 0, 0], ...    % Black edge for violins (RGB triplet)
    'BoxColor', [0, 0, 0], ...     % Black whisker and box lines (RGB triplet)
    'ShowData', true, ...          % Show data points
    'ShowMedian', true);           % Show median marker
title('Customized Violin Plot of GSI across Locations and SSP Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('GSI', 'FontSize', 12);
xlabel('Group (Location + SSP)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;



