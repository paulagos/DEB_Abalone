function get_plots_medium(simu, gen_number, legendLabels, colorIndex, foodQualities, SSPs, locations, h, i, j)
% Handles plotting when multiple locations and one SSP scenario OR multiple SSP scenarios and one location are selected

subPlotSeriesVec = foodQualities;
subPlotSeriesVar = h;

% Define line styles and colors for plotting

n_series = length(subPlotSeriesVec);  % Adjust the number of colors needed
colors = hsv(n_series);  % Create distinct colors using MATLAB's hsv function


if length(locations) > 1
    
    % Case 1: Multiple locations, one SSP scenario
    fprintf('Plotting for Location: %s, SSP: %s\n', locations{j}, SSPs{1});
    
    subPlotRowVec  = locations;
    subPlotRowVar  = j;
    
    figure(1);
    sgtitle(sprintf('Observables for SSP: %s across Multiple Locations', SSPs{1}));
    
    if ( i == 1 && j == 1 && h == 1 )
        % add a bit space to the figure. only do this once
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
    end
    
    for gen = 1:gen_number
        
        % Extract time and observable variables
        t = simu.tEVHR(:, 1, gen);  % Time vector
        L_w = simu.obs(:, 3, gen);  % Physical length
        Ww_V = simu.obs(:, 4, gen);  % Wet-weight
        GSI = simu.obs(:, 8, gen);  % Gonado Somatic Index (GSI)
        
        %% state variables
        E_H = simu.tEVHR(:,4, gen); % maturity
        E_R = simu.tEVHR(:,5, gen); % reproduction buffer
        
        %% Assign dates to simulation time span
        years = t(end)/365;
        
        % Get today and the start of the current year
        d = datetime(2019+(gen-1)*simu.tm/365,01,1);
        s = dateshift(d, 'start', 'year');
        t0Date = s + simu.tb - 1;
        tEndDate = s + years * 366 + simu.tb - 1; % using 366 instead of 365 as the date functions include leap years, and dont want to wind up with too few dates once we cut them out
        
        sim_dates = t0Date:tEndDate;
        sim_dates = sim_dates(~(sim_dates.Month == 2 & sim_dates.Day == 29));  % removing the leap days as we do not account for those in sim (all our years are 365 days)
        fprintf('Plotting for Location: %s, SSP: %s, cohort: %i\n', locations{j}, SSPs{i}, gen);
        
        %% find spawning dates indices
        i_sp = find(and((E_R == 0),(E_H >= simu.zoom^3 * simu.par.E_Hp)));
        
        %% Applyingthe dates to i_sp
        Spawning_dates = sim_dates(i_sp);
        
        %% ================ Create subplots of obervables for each location ================
        
        subplot(length(subPlotRowVec), 3, (subPlotRowVar-1)*3 + 1); hold on;
        plot(t/365, L_w, 'Linewidth', 1.5, 'LineStyle', '-', 'Color', colors(colorIndex,:) );
        title(sprintf('Length (Location: %s)', locations{j}));
        
        subplot(length(subPlotRowVec), 3, (subPlotRowVar-1)*3 + 2); hold on;
        plot(t/365, Ww_V, 'Linewidth', 1.5, 'LineStyle', '-', 'Color', colors(colorIndex,:) );
        title(sprintf('Weight (Location: %s)', locations{j}));
        
        subplot(length(subPlotRowVec), 3, (subPlotRowVar-1)*3 + 3); hold on;
        plot(t/365, GSI, 'Linewidth', 1.5, 'LineStyle', '-', 'Color', colors(colorIndex,:) );
        title(sprintf('GSI (Location: %s)', locations{j}));
        
        
        if ( gen == gen_number)
            
            lineChildren  = findobj(gca, 'Type', 'line');
            legendHandles = lineChildren(gen*subPlotSeriesVar:-gen_number:gen);
            
            % add legend
            Lgnd = legend(legendHandles ,legendLabels);
            Lgnd.Position(1) = 0.02;
            Lgnd.Position(2) = 0.5;
        end
        
    end  
        
        %% ===================== Fig.2 Spawning events  =====================
        
        figure(2)
        
        sgtitle('Spawning Events for Multiple Locations and SSP Scenarios');
        
        if ( i == 1 && j == 1 && h == 1 && gen == 1)
            % add a bit space to the figure. only do this once
            fig = gcf;
            fig.Position(3) = fig.Position(3) + 250;
        end
        
        spa_2 = Spawning_dates;
        spa_2.Year = 2020;
        
        jitterAmountDays = 1;  % Adjust this for more or less jitter on the dates
        jitterAmountYears = 0.5; % Adjust this for more or less jitter on the years
        
        datesJittered = spa_2 + (rand(size(spa_2)) - 0.5) * jitterAmountDays;
        yearsJittered = Spawning_dates.Year + (rand(size(Spawning_dates.Year)) - 0.5) * jitterAmountYears;
        
        plot(yearsJittered, datesJittered, 'o', 'Color', colors(colorIndex,:), 'MarkerFaceColor', colors(colorIndex,:))
        if (~isempty(Spawning_dates.Year))
            ytickformat('dd-MMM')
            tstart = datetime(2020,04,01);
            tend = datetime(2020,10,31);
            ylim([tstart tend])
        end
        title(sprintf('Location: %s, SSP: %s', locations{j}, SSPs{i}));
        ylabel('Dates')
        xlabel('Time (years)')
        
    
else
    
    
    % Case 2: One location, multiple SSP scenarios
    fprintf('Plotting for Location: %s, SSP: %s\n', locations{1}, SSPs{i});
    
    subPlotRowVec  = SSPs;
    subPlotRowVar  = i;
    
    figure(1);
    sgtitle(sprintf('Observables for Location: %s across Multiple SSP Scenarios', locations{1}));
    
    if ( i == 1 && j == 1 && h == 1 )
        % add a bit space to the figure. only do this once
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
    end
    
    for gen = 1:gen_number
        
        % Extract time and variables
        t = simu.tEVHR(:, 1, gen);  % Time vector
        L_w = simu.obs(:, 3, gen);  % Physical length
        Ww_V = simu.obs(:, 4, gen);  % Wet-weight
        GSI = simu.obs(:, 8, gen);  % Gonado Somatic Index (GSI)
        
        % ================ Create subplots of obervables for each location ================
        
        subplot(length(subPlotRowVec), 3, (subPlotRowVar-1)*3 + 1); hold on;
        plot(t/365, L_w, 'Linewidth', 1.5, 'LineStyle', '-', 'Color', colors(colorIndex,:) );
        title(sprintf('Length (SSP: %s)', SSPs{i}));
        
        
        subplot(length(subPlotRowVec), 3, (subPlotRowVar-1)*3 + 2); hold on;
        plot(t/365, Ww_V, 'Linewidth', 1.5, 'LineStyle', '-', 'Color', colors(colorIndex,:));
        title(sprintf('Weight (SSP: %s)', SSPs{i}));
        
        subplot(length(subPlotRowVec), 3, (subPlotRowVar-1)*3 + 3); hold on;
        plot(t/365, GSI, 'Linewidth', 1.5, 'LineStyle', '-', 'Color', colors(colorIndex,:) );
        title(sprintf('GSI (SSP: %s)', SSPs{i}));
        
        if ( gen == gen_number)
            
            lineChildren  = findobj(gca, 'Type', 'line');
            legendHandles = lineChildren(gen*subPlotSeriesVar:-gen_number:gen);
            
            % add legend
            Lgnd = legend(legendHandles ,legendLabels);
            Lgnd.Position(1) = 0.02;
            Lgnd.Position(2) = 0.5;
        end
        
        
    end
end

end
