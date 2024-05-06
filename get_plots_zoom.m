function get_plots_zoom(simu, SSPs, locations, gen_number, i, SSP, location)

if simu.zoom == 1
    colors = 'k';
    lineWidth = 1.6;
else
    colors = [.7 .7 .7];
    lineWidth = 1;
end

for gen = 1:gen_number
    
    %% observable quantities
    t    = simu.tEVHR(:,1,gen);  % Simulation time
    L_w  = simu.obs(:,3, gen);   % Physical length
    Ww_V = simu.obs(:,4, gen);   % Total wet-weight
    
    
    %%%% plot growth at with zoom factor %%%%%
    figure (6)
    subplot(SSPs, locations, i)
    hold on
    plot(t/365, Ww_V, 'Color', colors, 'LineWidth', lineWidth)
    title([SSP,' ', location])
    xlabel('Time (years)')
    ylabel('Total wet-weight (g)')
    
    figure (7)
    subplot(SSPs, locations, i)
    hold on
    plot(t/365, L_w, 'Color', colors, 'LineWidth', lineWidth) % cm, physical length
    title([SSP,' ', location])
    xlabel('Time (years)')
    ylabel('Total length (cm)')
    
end

end
    
 