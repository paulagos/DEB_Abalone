function get_plots(simu, h, i, gen_number, legendLabels, colorIndex, location, mode, SSP, totalLocations, totalSSP)
%---------------------------------------------------------------
% simu: simulation data structure
% h, i: indices for the selected options
% gen_number: number of generations being simulated
% legendLabels: dynamic labels for the legend (e.g., SSP_scenarios, Locations)
% colorIndex: index used to assign unique colors
% location: the current location being simulated
% mode: specifies the mode of plotting (e.g., default, locations_together, etc.)
% SSP: the current SSP scenario being processed
% totalLocations: total number of locations selected
% totalSSP: total number of SSP scenarios selected

% Display the running SSP scenario and location in the Command Window
fprintf('Running simulation for SSP: %s, Location: %s\n', SSP, location);

% Define dynamic line styles and color palette
line_styles = {'-', '--', ':', '-.'}; % Line styles 
n_colors = 20;  % Define how many distinct colors you need
colors = hsv(n_colors);  % Create distinct colors using MATLAB's hsv function

% Dynamically adjust figure size based on the number of locations
figureHeight = 300 + (150 * totalLocations);  % Base height with adjustment for more locations
figureWidth = 900;  % Fixed width to accommodate 3 columns (length, weight, GSI)

% Create Figure 1 for observables with dynamic size
figure(1);
set(gcf, 'Position', [100, 100, figureWidth, figureHeight]);  % Set dynamic size
sgtitle('Observables for Multiple Locations and SSP Scenarios');

% Loop through the generations
for gen = 1:gen_number
    % Extract the time vector `t` from the simulation data
    t = simu.tEVHR(:,1, gen);  % Assuming time is in the first column of `tEVHR`
    
    % Extract other variables for plotting
    L_w = simu.obs(:,3, gen);  % Physical length
    Ww_V = simu.obs(:,4, gen);  % Total wet-weight
    GSI = simu.obs(:,8, gen);  % Gonado Somatic Index (GSI)
    
    % ===================== Figure 1: Observables =====================
    % Create a 3-column figure: rows for locations, columns for observables
    % Rows = totalLocations, Columns = 3 (Length, Weight, GSI)
    
    % Row index for this location (h = location index)
    rowOffset = (h - 1) * 3;
    
    % Column 1: Physical Length (L_w)
    subplot(totalLocations, 3, rowOffset + 1); hold on;
    plot(t/365, L_w, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title(sprintf('Physical Length (Location: %s)', location));
    xlabel('Time (years)');
    ylabel('Length (cm)');

    % Column 2: Wet Weight (Ww_V)
    subplot(totalLocations, 3, rowOffset + 2); hold on;
    plot(t/365, Ww_V, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title(sprintf('Wet Weight (Location: %s)', location));
    xlabel('Time (years)');
    ylabel('Weight (g)');

    % Column 3: Gonado Somatic Index (GSI)
    subplot(totalLocations, 3, rowOffset + 3); hold on;
    plot(t/365, GSI, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title(sprintf('GSI (Location: %s)', location));
    xlabel('Time (years)');
    ylabel('GSI');

    legend(legendLabels, 'Location', 'best');

    % ===================== Separate Figures for Fluxes =====================
    % Create separate figures for each location and plot all SSPs within it
    figure('Name', ['Model Fluxes for Location: ', location], 'NumberTitle', 'off');
    sgtitle(['Model Fluxes for Location: ', location]);

    % Plot model fluxes for the current location and SSPs
    plot_fluxes(t, simu, gen, i, colorIndex, line_styles, colors, n_colors, legendLabels);

    % ===================== Separate Figures for Environmental Variables =====================
    figure('Name', ['Environmental Variables for Location: ', location], 'NumberTitle', 'off');
    sgtitle(['Environmental Variables for Location: ', location]);

    % Plot environmental variables for the current location and SSPs
    plot_environmental(t, simu, gen, i, colorIndex, line_styles, colors, n_colors, legendLabels);
end
end

% Helper function to plot fluxes
function plot_fluxes(t, simu, gen, i, colorIndex, line_styles, colors, n_colors, legendLabels)
    % Model fluxes (pX, pA, pM, etc.)
    pX = simu.tFluxes(:,2, gen);
    pA = simu.tFluxes(:,3, gen);
    pM = simu.tFluxes(:,4, gen);
    pJ = simu.tFluxes(:,5, gen);
    pC = simu.tFluxes(:,6, gen);
    pG = simu.tFluxes(:,7, gen);
    pR = simu.tFluxes(:,8, gen);
    pD = simu.tFluxes(:,9, gen);

    % Feeding flux (pX)
    subplot(3, 3, 1); hold on;
    plot(t, pX, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('Feeding (pX)');
    xlabel('Time (d)');
    ylabel('pX');
    
    % Assimilation flux (pA)
    subplot(3, 3, 2); hold on;
    plot(t, pA, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('Assimilation (pA)');
    xlabel('Time (d)');
    ylabel('pA');
    
    % Somatic Maintenance flux (pM)
    subplot(3, 3, 3); hold on;
    plot(t, pM, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('Somatic Maintenance (pM)');
    xlabel('Time (d)');
    ylabel('pM');
    
    % Maturity flux (pJ)
    subplot(3, 3, 4); hold on;
    plot(t, pJ, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('Maturity (pJ)');
    xlabel('Time (d)');
    ylabel('pJ');
    
    % Reserve Mobilization (pC)
    subplot(3, 3, 5); hold on;
    plot(t, pC, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('Reserve Mobilization (pC)');
    xlabel('Time (d)');
    ylabel('pC');
    
    % Growth (pG)
    subplot(3, 3, 6); hold on;
    plot(t, pG, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('Growth (pG)');
    xlabel('Time (d)');
    ylabel('pG');
    
    % Reproduction (pR)
    subplot(3, 3, 7); hold on;
    plot(t, pR, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('Reproduction (pR)');
    xlabel('Time (d)');
    ylabel('pR');
    
    % Dissipation (pD)
    subplot(3, 3, 8); hold on;
    plot(t, pD, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('Dissipation (pD)');
    xlabel('Time (d)');
    ylabel('pD');
    
    legend(legendLabels, 'Location', 'best');
end

% Helper function to plot environmental variables
function plot_environmental(t, simu, gen, i, colorIndex, line_styles, colors, n_colors, legendLabels)
    % Environmental variables
    X = simu.tEVHR(:,2, gen);  % Placeholder for food biomass (X)
    T = simu.tempData(t+1) - 273.15;  % Temperature in Celsius
    pCO2 = simu.pco2Data(t+1);  % pCO2 levels

    % Food Biomass (X)
    subplot(3, 1, 1); hold on;
    plot(t/365, X, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('Food Biomass (X)');
    xlabel('Time (years)');
    ylabel('Biomass');
    
    % Temperature (T)
    subplot(3, 1, 2); hold on;
    plot(t/365, T, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('Temperature (T)');
    xlabel('Time (years)');
    ylabel('Temperature (°C)');
    
    % pCO2
    subplot(3, 1, 3); hold on;
    plot(t/365, pCO2, 'Linewidth', 1.5, 'LineStyle', line_styles{mod(i-1, length(line_styles))+1}, 'Color', colors(mod(colorIndex, n_colors) + 1, :));
    title('pCO2');
    xlabel('Time (years)');
    ylabel('pCO2');
    
    legend(legendLabels, 'Location', 'best');
end


