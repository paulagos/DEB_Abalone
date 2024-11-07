function get_plots_simple(simu, gen_number, legendLabels, colorIndex, foodQualities, location, SSP, h)
    % Handles plotting when only one SSP scenario and one location are selected
    
    fprintf('Running simulation for SSP: %s, Location: %s\n', SSP, location);
    
    subPlotSeriesVec = foodQualities;
    subPlotSeriesVar = h;
    
    % Define line styles and colors for plotting
    line_styles = {'-', '--', ':', '-.'};
    n_series = length(subPlotSeriesVec);  % Adjust the number of colors needed
    colors   = hsv(n_series);             % Create distinct colors using MATLAB's hsv function


    for gen = 1:gen_number
        
        % Extract time and variables
        t = simu.tEVHR(:, 1, gen);  % Time vector
        L_w = simu.obs(:, 3, gen);  % Physical length
        Ww_V = simu.obs(:, 4, gen);  % Wet-weight
        GSI = simu.obs(:, 8, gen);  % Gonado Somatic Index (GSI)
        
        % Plot Observables in a single figure
        figure(1); 
        sgtitle(sprintf('Observables for SSP: %s, Location: %s', SSP, location));
        
        if ( h == 1 && gen == 1)
            % add a bit space to the figure. only do this once
            fig = gcf;
            fig.Position(3) = fig.Position(3) + 250;
        end

        subplot(1, 3, 1); hold on
        plot(t/365, L_w, 'Linewidth', 1.5, 'LineStyle', line_styles{1}, 'Color', colors(colorIndex, :));
        title('Physical Length');
        xlabel('Time (years)');
        ylabel('Length (cm)');

        subplot(1, 3, 2); hold on
        plot(t/365, Ww_V, 'Linewidth', 1.5, 'LineStyle', line_styles{1}, 'Color', colors(colorIndex, :));
        title('Wet Weight');
        xlabel('Time (years)');
        ylabel('Weight (g)');

        subplot(1, 3, 3); hold on
        plot(t/365, GSI, 'Linewidth', 1.5, 'LineStyle', line_styles{1}, 'Color', colors(colorIndex, :));
        title('GSI');
        xlabel('Time (years)');
        ylabel('GSI');
        
        if ( gen == gen_number)
   
            lineChildren  = findobj(gca, 'Type', 'line');
            legendHandles = lineChildren(gen*subPlotSeriesVar:-gen_number:gen);
            
            % add legend
            Lgnd = legend(legendHandles ,legendLabels);
            Lgnd.Position(1) = 0.02;
            Lgnd.Position(2) = 0.5;
        end

        % Similar approach can be used to plot fluxes and environmental variables
    end
    
end