function get_plots_complex(simu, gen_number, legendLabels, colorIndex, foodQualities, SSPs, locations, zooms, h, i, j, m)


%% Handles plotting when multiple locations and multiple SSP scenarios are selected
subPlotRowVec  = SSPs;
subPlotRowVar  = i;

subPlotColVec  = locations;
subPlotColVar  = j;

subPlotSeriesVec = foodQualities;
subPlotSeriesVar = h;

currentSubplot = (subPlotRowVar-1)*length(subPlotColVec) + subPlotColVar;


%% Define line styles and colors for plotting

n_series   = length(subPlotSeriesVec);  % Adjust the number of colors needed
colors     = hsv(n_series);  % Create distinct colors using MATLAB's hsv function
colors_env = hsv(length(SSPs));

line_styles = {'-', '--', ':', '-.'}; % Line styles
myLineStyle = line_styles{h};

if zooms(m) == 1
    myLineWidth = 1.5;
    myColor = colors(h,:);
    myColorEnv = colors_env(i,:);
else
    myLineWidth = 0.5;
    myColor = [.7 .7 .7];
    myColorEnv = [.7 .7 .7];
end

%% Define data to plot in time

for gen = 1:gen_number

 %% environmental forcing
    t = simu.tEVHR(:,1, m, gen);
    
    if simu.envX == 1
        X = repmat(simu.Xinit, length(t),1);
    else
        X = food(t,simu);
    end
    
    if simu.envT == 1
        T = repmat(simu.Tinit, length(t), 1) - 273.15; % in DegC
    else
        T = simu.tempData(t+1);
    end
    
    if simu.envpCO2 == 1
        pCO2 = repmat(simu.pCO2init, length(t), 1);
    else
        pCO2 = simu.pco2Data(t+1);
    end
    
    %% state variables
    E   = simu.tEVHR(:,2, m, gen); % Energy reserves
    V   = simu.tEVHR(:,3, m, gen); % structural volume
    E_H = simu.tEVHR(:,4, m, gen); % maturity
    E_R = simu.tEVHR(:,5, m, gen); % reproduction buffer
    
    %% Model fluxes
    pX  = simu.tFluxes(:,2, m, gen);
    pA  = simu.tFluxes(:,3, m, gen);
    pM  = simu.tFluxes(:,4, m, gen);
    pJ  = simu.tFluxes(:,5, m, gen);
    pC  = simu.tFluxes(:,6, m, gen);
    pG  = simu.tFluxes(:,7, m, gen);
    pR  = simu.tFluxes(:,8, m, gen);
    pD  = simu.tFluxes(:,9, m, gen);
    
    %% Observables variables
    L_w  = simu.obs(:, 3, m, gen);    % Physical length
    Ww_V = simu.obs(:, 4, m, gen);    % Wet-weight
    GSI  = simu.obs(:, 8, m, gen);    % GSI
    
        %% Functional response
    f  = X ./ (X + simu.K); % Problem with K, SHOULD NOT COME FROM THE ESTIMATION
    
    %% Assign dates to simulation time span
    years = t(end)/365;
    
    % Get today and the start of the current year
    d = datetime(2019+(gen-1)*simu.tm/365,01,1);
    s = dateshift(d, 'start', 'year');
    t0Date = s + simu.tb - 1;
    tEndDate = s + years * 366 + simu.tb - 1; % using 366 instead of 365 as the date functions include leap years, and dont want to wind up with too few dates once we cut them out
    
    sim_dates = t0Date:tEndDate;
    sim_dates = sim_dates(~(sim_dates.Month == 2 & sim_dates.Day == 29));  % removing the leap days as we do not account for those in sim (all our years are 365 days)
    
    fprintf('Plotting for Food Quality: %.2f, Location: %s, SSP: %s, Zoom factor:, %.2f, cohort: %i\n', foodQualities(h), locations{j}, SSPs{i}, zooms(m), gen);
    
     %% find spawning dates indices
    i_sp = and((E_R == 0),(E_H >= simu.zoom^3 * simu.par.E_Hp));
    
     %% Applyingthe dates to i_sp
    Spawning_dates = sim_dates(i_sp);
    
    %% Create subplots for each location and SSP combination
 
    % ===================== Fig.1 Total lenghts  =====================
    figure(1);
    
    sgtitle(sprintf('Total Length for Multiple Locations and SSP Scenarios'));
    
    if ( i == 1 && j == 1 && h == 1 && m == 1 && gen == 1)
         % add a bit space to the figure. only do this once
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
    end
    
    subplot(length(subPlotRowVec), length(subPlotColVec), currentSubplot); hold on;
    
    plot(t/365, L_w, 'Linewidth', myLineWidth, 'LineStyle', '-', 'Color', myColor );
    title(sprintf('Location: %s, SSP: %s', locations{j}, SSPs{i}));
    ylabel('Total length (cm)')
    xlabel('Time (years)')
    
    
    if ( currentSubplot == 1 && m == length(zooms) && gen == gen_number)
       
        lineChildren  = findobj(gca, 'Type', 'line');
        %legendHandles = lineChildren(gen*subPlotSeriesVar:-gen_number:gen);
        legendHandles = lineChildren(end:-1:1);
        legendHandles = legendHandles(1:length(zooms)*gen_number:h*length(zooms)*gen_number);
      
        % add legend
        Lgnd = legend(legendHandles ,legendLabels);
        Lgnd.Position(1) = 0.02;
        Lgnd.Position(2) = 0.5;
    end
    
    % ===================== Fig.2 Shell free wet-weights  =====================
    figure(2);
    
    sgtitle('Wet Weight for Multiple Locations and SSP Scenarios');
    
    if ( i == 1 && j == 1 && h == 1 && m == 1 && gen == 1)
        % add a bit space to the figure. only do this once
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
    end
    
    subplot(length(subPlotRowVec), length(subPlotColVec), currentSubplot); hold on;
    
    plot(t/365, Ww_V, 'Linewidth', myLineWidth, 'LineStyle', '-', 'Color', myColor );
    title(sprintf('Location: %s, SSP: %s', locations{j}, SSPs{i}));
    ylabel('Shell free wet weight (g)')
    xlabel('Time (years)')
    
    if ( currentSubplot == 1 && m == length(zooms) && gen == gen_number)
       
        lineChildren  = findobj(gca, 'Type', 'line');
        %legendHandles = lineChildren(gen*subPlotSeriesVar:-gen_number:gen);
        legendHandles = lineChildren(end:-1:1);
        legendHandles = legendHandles(1:length(zooms)*gen_number:h*length(zooms)*gen_number);
      
        % add legend
        Lgnd = legend(legendHandles ,legendLabels);
        Lgnd.Position(1) = 0.02;
        Lgnd.Position(2) = 0.5;
    end
    
    % ===================== Fig.3 GSI  =====================
    figure(3);
    
    sgtitle('GSI for Multiple Locations and SSP Scenarios');
    
    if ( i == 1 && j == 1 && h == 1 && m == 1 && gen == 1)
        % add a bit space to the figure. only do this once
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
    end
    
    subplot(length(SSPs), length(locations), currentSubplot); hold on;
    
    plot(t/365, GSI, 'Linewidth', myLineWidth, 'LineStyle', '-', 'Color', myColor );
    title(sprintf('Location: %s, SSP: %s', locations{j}, SSPs{i}));
    ylabel('Gonado somatic index (GSI)')
    xlabel('Time (years)')
    
    if ( currentSubplot == 1 && m == length(zooms) && gen == gen_number)
       
        lineChildren  = findobj(gca, 'Type', 'line');
        %legendHandles = lineChildren(gen*subPlotSeriesVar:-gen_number:gen);
        legendHandles = lineChildren(end:-1:1);
        legendHandles = legendHandles(1:length(zooms)*gen_number:h*length(zooms)*gen_number);
      
        % add legend
        Lgnd = legend(legendHandles ,legendLabels);
        Lgnd.Position(1) = 0.02;
        Lgnd.Position(2) = 0.5;
    end
    
    
    % ===================== Fig.4 Spawning events  =====================
    figure(4)
    
    sgtitle('Spawning Events for Multiple Locations and SSP Scenarios');
    
    if ( i == 1 && j == 1 && h == 1 && m == 1 && gen == 1)
        % add a bit space to the figure. only do this once
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
    end

    subplot(length(SSPs), length(locations), currentSubplot); hold on;
    
    spa_2 = Spawning_dates;
    spa_2.Year = 2020;
    
    jitterAmountDays = 1;  % Adjust this for more or less jitter on the dates
    jitterAmountYears = 0.5; % Adjust this for more or less jitter on the years

    datesJittered = spa_2 + (rand(size(spa_2)) - 0.5) * jitterAmountDays;
    yearsJittered = Spawning_dates.Year + (rand(size(Spawning_dates.Year)) - 0.5) * jitterAmountYears;
    
    plot(yearsJittered, datesJittered, 'o', 'Color', colors(h,:), 'MarkerFaceColor', colors(h,:))
    if (~isempty(Spawning_dates.Year))
        ytickformat('dd-MMM')
        tstart = datetime(2020,04,01);
        tend = datetime(2020,10,31);
        ylim([tstart tend])
    end
    title(sprintf('Location: %s, SSP: %s', locations{j}, SSPs{i}));
    ylabel('Dates')
    xlabel('Time (years)')
    
    
    %%% Find lables for each food quality %%%%
    
    if currentSubplot == 1 && m == length(zooms) && gen == gen_number
        lineChildren = findobj(gca, 'Type', 'line');
        
        % Check if lineChildren is non-empty before creating a legend
        if ~isempty(lineChildren)
            % Reverse the order to match the original intent
            legendHandles = lineChildren(end:-1:1);
            
            % Indexing with safety to avoid "Index exceeds the number of elements" errors
            legendHandles = legendHandles(1:min(length(legendHandles), length(zooms)*gen_number));
            
            % Add legend if there are valid handles
            Lgnd = legend(legendHandles, legendLabels);
            Lgnd.Position(1) = 0.02;
            Lgnd.Position(2) = 0.5;
        else
            warning('No line objects found for legend creation in subplot.');
        end
    end
    
%     if ( currentSubplot == 1 && m == length(zooms) && gen == gen_number)
%        
%         lineChildren  = findobj(gca, 'Type', 'line');
%         legendHandles = lineChildren(gen*subPlotSeriesVar:-gen_number:gen);
%         legendHandles = lineChildren(end:-1:1);
%         legendHandles = legendHandles(1:length(zooms)*gen_number:h*length(zooms)*gen_number);
%       
%         % add legend
%         Lgnd = legend(legendHandles ,legendLabels);
%         Lgnd.Position(1) = 0.02;
%         Lgnd.Position(2) = 0.5;
%     end
%     
    
    % ===================== Fig 5: Fluxes =====================
    figure(5);
    
    sgtitle('DEB fluxes for Multiple Locations and SSP Scenarios');
    
    if ( i == 1 && j == 1 && h == 1 && m == 1 && gen == 1)
        % add a bit space to the figure. only do this once
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
    end
    
    subplot(8, length(locations), j); hold on;
    %subplot(3, 3, 1) % Feeding (pX)
    hold on;
    plot(t, pX, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title(sprintf(' Feeding (pX)  Location: %s', locations{j} ));
    xlabel('Time (d)');
    ylabel('pX');
    
    subplot(8, length(locations), length(locations)+j); hold on;
    %subplot(3, 3, 2) % Assimilation (pA)
    hold on;
    plot(t, pA, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv);
    title('Assimilation (pA)');
    xlabel('Time (d)');
    ylabel('pA');
    
    subplot(8, length(locations), 2*length(locations)+j); hold on;
    %subplot(3, 3, 3) % Somatic Maintenance (pM)
    hold on;
    plot(t, pM, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title('Somatic Maintenance (pM)');
    xlabel('Time (d)');
    ylabel('pM');
    
    subplot(8, length(locations), 3*length(locations)+j); hold on;
    %subplot(3, 3, 4) % Maturity Maintenance (pJ)
    hold on;
    plot(t, pJ, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title('Maturity Maintenance (pJ)');
    xlabel('Time (d)');
    ylabel('pJ');
    
    subplot(8, length(locations), 4*length(locations)+j); hold on;
    %subplot(3, 3, 5) % Reserve Mobilization (pC)
    hold on;
    plot(t, pC, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title('Reserve Mobilization (pC)');
    xlabel('Time (d)');
    ylabel('pC');
    
    subplot(8, length(locations), 5*length(locations)+j); hold on;
    %subplot(3, 3, 6) % Growth (pG)
    hold on;
    plot(t, pG, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title('Growth (pG)');
    xlabel('Time (d)');
    ylabel('pG');
    
    subplot(8, length(locations), 6*length(locations)+j); hold on;
    %subplot(3, 3, 7) % Reproduction (pR)
    hold on;
    plot(t, pR, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title('Reproduction (pR)');
    xlabel('Time (d)');
    ylabel('pR');
    
    subplot(8, length(locations), 7*length(locations)+j); hold on;
    %subplot(3, 3, 8) % Dissipation (pD)
    hold on;
    plot(t, pD, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title('Dissipation (pD)');
    xlabel('Time (d)');
    ylabel('pD');
    legend()
    
    
    % ===================== Fig.6 State variables  =====================
    figure(6);
    
    sgtitle('DEB State variables for Multiple Locations and SSP Scenarios');
    
    if ( i == 1 && j == 1 && h == 1 && m == 1 && gen == 1)
        % add a bit space to the figure. only do this once
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
    end
    
    subplot(5, length(locations), j); hold on;
    hold on;
    plot(t, E, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title(sprintf('E  Location: %s', locations{j} ));
    xlabel('Time (d)');
    ylabel('E (J)');
    
    subplot(5, length(locations), length(locations)+j); hold on;
    hold on;
    plot(t, V, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title(sprintf('V  Location: %s', locations{j} ));
    xlabel('Time (d)');
    ylabel('V (cm3)');
    
    subplot(5, length(locations), 2*length(locations)+j); hold on;
    hold on;
    plot(t, E_R, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title(sprintf('E_R  Location: %s', locations{j} ));
    xlabel('Time (d)');
    ylabel('E_R (J)');
    
    subplot(5, length(locations), 3*length(locations)+j); hold on;
    hold on;
    plot(t, E_H, 'Linewidth', myLineWidth, 'LineStyle', myLineStyle, 'Color', myColorEnv );
    title(sprintf('E_H  Location: %s', locations{j} ));
    xlabel('Time (d)');
    ylabel('E_H (J)');
    

    % ===================== Fig.7 distribution of spawning events  =====================
    figure(7)
    
    sgtitle('Frequency distribution for Multiple Locations and SSP Scenarios');
    
    if ( i == 1 && j == 1 && h == 1 && m == 1 && gen == 1)
        % add a bit space to the figure. only do this once
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
    end
    
    subplot(length(SSPs), length(locations), currentSubplot); hold on;
    
    spa_2 = Spawning_dates;  % Use the Spawning_dates for histogram
    
    % Convert Spawning_dates to a numeric format (e.g., days since start of the year)
    numeric_spawning_dates = datenum(spa_2);
    
    % Plot histogram without 'FaceColor' (as it's incompatible with 'stairs' DisplayStyle)
    num_bins = 5;  % Set to a higher number to make the columns narrower
    histogram(numeric_spawning_dates, num_bins,'FaceColor', 'blue', 'Normalization', 'count', 'FaceAlpha', 0.5);
    
    % Add labels and title
    ylabel('Frequency');
    xlabel('Spawning Dates');
    datetick('x', 'dd-mm-yy', 'keepticks');  % Format x-axis as dates
    title(sprintf('Location: %s, SSP: %s', locations{j}, SSPs{i}));
    
    
    % ===================== Fig. 8: Environmental Variables =====================
    figure (8)
    
    sgtitle('Environmental Variables'); % Overall title for environmental variables
    
    if ( i == 1 && j == 1 && h == 1 && m == 1 && gen == 1)
        % add a bit space to the figure. only do this once
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
    end
    
    if ( h == 1 && i == 1 && j == 1 && m == 1)
        subplot(3, length(subPlotColVec), subPlotColVar); hold on;
        plot(t/365, X, 'Linewidth', 1.5, 'LineStyle', myLineStyle, 'Color', myColor );
        title(sprintf('Food Biomass (X)'));
        xlabel('Time (years)');
        ylabel('Biomass');
        
        subplot(3, length(subPlotColVec), subPlotColVar+1); hold on;
        plot(t/365, f, 'Linewidth', 1.5, 'LineStyle', '-', 'Color', colors_env(i,:));
        title('Functional response');
        xlabel('Time (years)');
        ylabel('Functional response (f)');
    end
    
    if ( h == 1 && m == 1 )
        subplot(3, length(subPlotColVec), subPlotColVar + length(subPlotColVec)); hold on;
        plot(t/365, T, 'Linewidth', 1.5, 'LineStyle', '-', 'Color', colors_env(i,:));
        title('Temperature (T)');
        xlabel('Time (years)');
        ylabel('Temperature (°C)');
        
        subplot(3, length(subPlotColVec), subPlotColVar + 2 * length(subPlotColVec)); hold on;
        plot(t/365, pCO2, 'Linewidth', 1.5, 'LineStyle', '-', 'Color', colors_env(i,:));
        title('pCO2');
        xlabel('Time (years)');
        ylabel('pCO2');
        
        if (gen == gen_number && subPlotColVar == 1)
            
            lineChildren  = findobj(gca, 'Type', 'line');
            legendHandles = lineChildren(gen*subPlotRowVar:-gen_number:gen);
            
            legendLabelsEnv = subPlotRowVec(1:length(legendHandles)); 
            
            % add legend
            Lgnd = legend(legendHandles, legendLabelsEnv);
            Lgnd.Position(1) = 0.02;
            Lgnd.Position(2) = 0.5;
        end
        
    end
    
 

    
end

end