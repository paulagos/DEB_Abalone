%---------------------------------------------------------------
% Compute model predictions for one individual,
%   with standard DEB equations (Kooijman 2010)
%
% Parameter estimates for Haliotis tuberculata
%
% from the egg to adult stage
% one reproduction event every year
% (shrinking of structural volume is possible in starvation condition)
%
% simu : structure with individual features and model predictions
%
% calls :   init.m
%           indiv.m
%           get_fluxes.m
%           get_obs.m
%           get_plots.m
%
% 2013/03/15 - Laure Pecquerie
% 2023/03/12 - Paulo F. Lagos
%---------------------------------------------------------------

clear all
close all

RCP_scenarios = {'RCP1','RCP2','RCP3','RCP4'};
Food_quality  = [-0.4, -0.2, 0];
variable_env  = [1,2];
gen_number    = [2];

for h = 1:length(Food_quality)
    envFoodQ = Food_quality(h);
    
    for i = 1:length(RCP_scenarios)
        RCP = RCP_scenarios{i};
        
        for j = 1:length(variable_env)
            varEnv = variable_env(j);
            
            %% 1 - initialize time, parameters, etc
            simu = init(RCP, varEnv, envFoodQ, gen_number);
            
            for k = 1:gen_number
                
                %% 2 - calculate state variables
                simu.tEVHR(:,:,k) = indiv(simu, k);
                
                %% 3 - calculate model fluxes
                simu.tFluxes(:,:,k) = get_fluxes(simu, k);
                
                %% 4 - calculate observable quantities
                simu.obs(:,:,k) = get_obs(simu, k);
                
            end
            
            %% 5 - save data structure
            eval( ['simu_foodQ', num2str(h), '_' ,RCP, '_env', num2str(varEnv), '= simu;'] )
        
        end
    end
end


for h = 1:length(Food_quality)
    
    for i = 1:length(RCP_scenarios)
        
        RCP = RCP_scenarios{i};
        
        varEnv = variable_env(1); % for ploting  constant (1) or variable (2) environmental conditions
        %h = 1;
        
        eval( ['simu = simu_foodQ', num2str(h), '_' ,RCP, '_env', num2str(varEnv), ';'] )
        
        %% 6 - make plots
        get_plots(simu, i, h, gen_number);
        
    end
end