function get_plots(simu, j, h, gen_number)

% choose plots
doLarvalPlot = 0; % chose if do the larval plots 0 = false; 1 = true

% Construct the filename based on the scenario
filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';

% Load parameter values
res = load(filename);
par  = res.par;          % get primary parameters
cPar = parscomp_st(par); % get compound parameters

for gen = 1:gen_number
    
    % time , environmental forcing
    t = simu.tEVHR(:,1,gen);
    
    if simu.envX == 1
        X = repmat(simu.Xinit, length(t),1);
    else
        X = food(t);
    end
    
    if simu.envT == 1
        T = repmat(simu.Tinit, length(t), 1) - 273; % in DegC
    else
        T = temp(t) - 273.15;
    end
    
    if simu.envpCO2 == 1
        pCO2 = repmat(simu.pCO2init, length(t), 1);
    else
        pCO2 = acid(t);
    end
    
    % For period 2031 - 2050 scenarios (RCP2, RCP3 RCP4)
    if simu.RCP == "RCP2"  % period 2031 - 2050
        T = T + t * 0.64 / (365 * 10);
        pCO2 = pCO2 + t * 10 / (365 * 10);
    end
    
    if simu.RCP == "RCP3" % period 2031 - 2050
        T = T + t * 0.95 / (365 * 10);
        pCO2 = pCO2 + t * 15 / (365 * 10);
    end
    
    if simu.RCP == "RCP4" % period 2031 - 2050
        T = T + t * 1.5 / (365 * 10);
        pCO2 = pCO2 + t * 120 / (365 * 10);
    end
    
    
    %% state variables
    E   = simu.tEVHR(:,2, gen); % Energy reserves
    V   = simu.tEVHR(:,3, gen); % structural volume
    E_H = simu.tEVHR(:,4, gen); % maturity
    E_R = simu.tEVHR(:,5, gen); % reproduction buffer
    
    %% spawning date indices
    i_sp = find(and((E_R == 0),(E_H >= par.E_Hp)));
    i_sp = i_sp - 1 ; % look at the preceding line with E_R value before spawning
    
    %% temperature correction
    c_T2  = exp(par.T_A / par.T_ref - par.T_A/ 297.8);
    
%     s_A = exp(par.T_A/ par.T_ref - par.T_A ./ T);  % Arrhenius factor
%     
%     s_L_ratio = (1 + exp(par.T_AL/ par.T_ref - par.T_AL/ par.T_L)) ./ ...
%         (1 + exp(par.T_AL ./ T   - par.T_AL / par.T_L));
%     
%     s_H_ratio = (1 + exp(par.T_AH / par.T_H - par.T_AH/ par.T_ref)) ./ ...
%         (1 + exp(par.T_AH / par.T_H - par.T_AH ./ T  ));
%     
%     c_T = s_A .* ((T <= par.T_ref) .* s_L_ratio + (T > par.T_ref) .* s_H_ratio);
    
    %% Functional response
    f  = X ./ (X + cPar.K);
    
    %% Model fluxes
    pX  = simu.tFluxes(:,2, gen);
    pA  = simu.tFluxes(:,3, gen);
    pM  = simu.tFluxes(:,4, gen);
    pJ  = simu.tFluxes(:,5, gen);
    pC  = simu.tFluxes(:,6, gen);
    pG  = simu.tFluxes(:,7, gen);
    pR     = simu.tFluxes(:,8, gen);
    pD     = simu.tFluxes(:,9, gen);
    pLR    = simu.tFluxes(:,10, gen);
    pLV    = simu.tFluxes(:,11, gen);
    kap_pC = simu.tFluxes(:,12, gen);
    
    %% Larval age at birth and metamorphosis
    if doLarvalPlot
        pars_tj = [cPar.g; cPar.k; cPar.l_T; cPar.v_Hb; cPar.v_Hj; cPar.v_Hp];
        [t_j, ~, t_b, l_j, ~, l_b, ~, ~, ~, ~] = get_tj(pars_tj, f);
        
        aT_b = t_b / cPar.k_M / c_T2;    % d, age at birth at f and T
        L_b = cPar.L_m * l_b;             % cm, structural length at birth at f (post-metamorphosis)
        Lw_b = L_b / par.del_L;          % cm, total length at birth at f (post-metamorphosis)
        Ww_b = L_b^3 * (1 + cPar.w * cPar.g);  % g, wet-weight at birth
        
        tT_j = t_j / cPar.k_M / c_T2;    % d, age at metamorphosis
        L_j = cPar.L_m * l_j;             % cm, structural length at metam
        Lw_j = L_j / par.del_L;          % cm, total lenght at metamorphosis at f
        Ww_j = L_j^3 * (1 + 0.1 * cPar.w);  % g, wet-weight at metam
        
    end
    %% observable quantities
    e    = simu.obs(:,2, gen); % Scaled reserve energy
    L_w  = simu.obs(:,3, gen); % Physical length
    Ww_V = simu.obs(:,4, gen); % Total wet-weight
    E_w  = simu.obs(:,5, gen); % Total energy content
    F    = simu.obs(:,6, gen); % Fecundity
    G_Ww = simu.obs(:,7, gen); % Gonad wet weight
    
        %% graph 1
    Enfig = 1; % number of figures
    for i = 1:nfig
        figure(i)
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','in');
        set(gca,'Fontsize', 12)
        % Adjust the figure size, if needed
        set(gcf, 'Position', [100, 100, 800, 600]);
    end
    
    line_styles = {'-', '--', ':', '-.'}; % Line styles for IPCC scenarios
    colors = ["#77AC30", "#0072BD", "#7E2F8E" ,"#A2142F"];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Plots for model fluxes %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    subplot(3, 4, 1)
    hold on
    plot(t, pX, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Feeding')
    %legend('Current','RCP2.6','RCP4.5','RCP8.5', 'location', 'best'),
    xlabel('Time (d)')
    ylabel('pX')
    
    subplot(3, 4, 2)
    hold on
    plot(t, pA, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Assimilation')
    xlabel('Time (d)')
    ylabel('pA')
    
    subplot(3, 4, 3)
    hold on
    plot(t, pM, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Somatic maintenance')
    xlabel('Time (d)')
    ylabel('pM')
    
    subplot(3, 4, 4)
    hold on
    plot(t, pJ, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Maturity maintenance')
    xlabel('Time (d)')
    ylabel('pJ')
    
    subplot(3, 4, 5)
    hold on
    plot(t, pC, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Reserve mobilization')
    xlabel('Time (d)')
    ylabel('pC')
    
    subplot(3, 4, 6)
    hold on
    plot(t, pG, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Growth')
    xlabel('Time (d)')
    ylabel('pG')
    
    subplot(3, 4, 7)
    hold on
    plot(t, pR, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Reproduction')
    xlabel('Time (d)')
    ylabel('pR')
    
    subplot(3, 4, 8)
    hold on
    plot( t, pD, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Dissipation')
    xlabel('Time (d)')
    ylabel('pD')
    
    subplot(3, 4, 9)
    hold on
    plot( t, pLR , simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('pLR')
    xlabel('Time (d)')
    ylabel('pLR ')
    
    subplot(3, 4, 10)
    hold on
    plot( t, pLV, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('pLV')
    xlabel('Time (d)')
    ylabel('pLV')
    
    subplot(3, 4, 11)
    hold on
    plot( t, kap_pC, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title(' kap * pC')
    xlabel('Time (d)')
    ylabel('kap * pC')
    
    subplot(3, 4, 12)
    p1 = plot(NaN,NaN, 'Linewidth',1.5, 'Color', colors{1});
    hold on;
    p2 = plot(NaN,NaN, 'Linewidth',1.5, 'Color', colors{2});
    p3 = plot(NaN,NaN, 'Linewidth',1.5, 'Color', colors{3});
    p4 = plot(NaN,NaN, 'Linewidth',1.5, 'Color', colors{4});
    
    legend('Current','RCP2.6','RCP4.5','RCP8.5', 'location', 'best');

    set(gca, 'visible', 'off');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Plots for Enviromental conditions %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    subplot(3, 1, 1)
    hold on
    plot(t, T, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Temperature')
    xlabel('Time (years)')
    ylabel('T (°C)')
    axis([0 t(end) 0 25])
    
    subplot(3, 1, 2)
    hold on
    plot(t/365, f , simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Food')
    xlabel('Time (years)')
    ylabel('Food density')
    
    subplot(3, 1, 3)
    hold on
    plot(t/365, pCO2, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('pCO2 levels')
    xlabel('Time (years)')
    ylabel('pCO2')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Plots for State variables  %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(3)
    subplot(2, 4, 1)
    hold on
    plot(t/365, e, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})
    title('Scaled reserve density (e)')
    xlabel('Time (d)')
    ylabel('e (-)')
    
    subplot(2, 4, 2)
    hold on
    plot(t/365, E,  simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j}) % Joules, Reserve energy
    title('Reserve energy [E]')
    xlabel('Time (d)')
    ylabel('E (Joules)')
    
    subplot(2, 4, 3)
    hold on
    plot(t/365, V, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j}) % cm^3, Structural volume
    title('Structural volume')
    xlabel('Time (d)')
    ylabel('V (cm^3)')
    
    subplot(2, 4, 4)
    hold on
    plot(t/365, E_H, simu.col, 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})  % J, cumulated energy inversted into maturity
    title('Energy invested into maturation')
    xlabel('Time (d)')
    ylabel('E_H (J)')
    
    subplot(2, 4, 5)
    hold on
    plot(t/365, E_R, simu.col , 'Linewidth',1.5, 'LineStyle', line_styles{h}, 'Color', colors{j}) % Joules, Reproduction buffer
    title('Reproduction buffer')
    xlabel('Time (d)')
    ylabel('E_R (J)')
    
    subplot(2, 4, 6)
    hold on
    plot(t/365, E_w, simu.col , 'Linewidth', 1.5, 'LineStyle', line_styles{h}, 'Color', colors{j}) % J/g, Total energy content
    title('Energy content (J/g)')
    xlabel('Time (d)')
    ylabel('E_w (J.g^{-1})')
    
    saveas(gcf, 'fig_DEB_standard.jpg', 'jpg')
    % other ways to save the graph
    % saveas(gcf, 'fig_DEB_standard.eps', 'psc2')
    % print -djpeg90 fig_DEB_standard.jpg
    
    
    % graphs 2
    nfig = 1; % number of figures
    for i = 1:nfig
        figure(i)
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','in');
        set(gcf,'PaperPosition',[0 0 14 10]);
        set(gca,'Fontsize', 14)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Plots for observable quantities %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%% Total lenghts %%%%%%%%
    figure(4)
    subplot(1, 4, 1)
    hold on
    plot(t/365, L_w, simu.col, 'Linewidth', 1.5, 'LineStyle', line_styles{h}, 'Color', colors{j}) % cm, physical length
    xlabel('Time (years)')
    ylabel('Total length (cm)')
    v = axis;
    text((v(2)-v(1)) * 0.07 + v(1), (v(4)-v(3)) *.9 + v(3), '(c)','Fontsize',12)
    
    %%%%%%%%% Fecundity %%%%%%%%
    subplot(1, 4, 2)
    hold on
    plot(t/365, F, simu.col, 'Linewidth', 1.5, 'LineStyle', line_styles{h}, 'Color', colors{j}) % Fecundity = egg number = kap_R * E_R / E_0
    xlabel('Time (years)')
    ylabel('Fecundity (# Eggs)')
    v = axis;
    text((v(2)-v(1)) * 0.07 + v(1), (v(4)-v(3)) *.9 + v(3), '(c)','Fontsize',12)
    
    %%%%%%%%% wet-weights %%%%%%%%
    subplot(1, 4, 3)
    hold on
    plot(t/365, Ww_V, simu.col, 'Linewidth', 1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})  % g, total wet-weight, same water content in structure and reserves
    xlabel('Time (years)')
    ylabel('Total wet-weight (g)')
    
    %%%%%% Gonad Weigth %%%%%%%%%
    subplot(1, 4, 4)
    hold on
    plot(t/365, G_Ww, simu.col, 'Linewidth', 1.5, 'LineStyle', line_styles{h}, 'Color', colors{j})  % g, total wet-weight, same water content in structure and reserves
    xlabel('Time (years)')
    ylabel('Gonad wet-weight (g)')
        
    
    if doLarvalPlot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Plots for larva %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Total lenghts %%%%%%%%
        figure(5)
        subplot(1, 2, 1)
        hold on
        plot(t, L_w, simu.col, 'Linewidth', 1.6) % cm, physical length
        xline (tT_j, '--.r', 'Predicted age at metamorphosis')
        axis([0 5 0 0.06])
        xlabel('Time (d)')
        ylabel('Larva length (cm)')
        
        %%%%%%%%% wet-weights %%%%%%%%
        subplot(1, 2, 2)
        hold on
        plot(t, Ww_V * 1000, simu.col, 'Linewidth', 1.6)  % g, total wet-weight, same water content in structure and reserves
        xline (tT_j, '--.r', 'Predicted age at metamorphosis')
        axis([0 5 0 0.004])
        xlabel('Time (d)')
        ylabel('Larva wet-weight (mg)')
        sgtitle('Larval growth')
    end
    
   
end
end

