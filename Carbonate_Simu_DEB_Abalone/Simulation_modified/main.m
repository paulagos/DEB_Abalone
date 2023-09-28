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
%           get_environmental_variables.m
%           get_obs.m
%           get_plots.m
%
% 2013/03/15 - Laure Pecquerie
% 2023/03/12 - Paulo F. Lagos
%---------------------------------------------------------------

clear all
%close all

RCP_scenarios = {'RCP1','RCP2','RCP3','RCP4'};
variable_env  = [1,2];

for i = 1:length(RCP_scenarios)
    RCP = RCP_scenarios{i};
    
    for j = 1:length(variable_env)
        varEnv = variable_env(j);
        
        %% 1 - initialize time, parameters, etc
        simu = init(RCP, varEnv);
        
        %% 2 - calculate state variables
        simu.tEVHRSO = indiv(simu);
        
        %% 3 - calculate model fluxes
        simu.tFluxes = get_fluxes(simu);
        
        %% 4 - calculate observable quantities
        simu.obs = get_obs(simu);
        
        %% 5 - save data structure
        eval(['simu_', RCP, '_env', num2str(varEnv), '= simu;'])
        
        
    end
end


for i = 1:length(RCP_scenarios)
    RCP = RCP_scenarios{i};
    
    varEnv = variable_env(1); % for ploting  constant (1) or variable (2) conditions
    eval(['simu = simu_', RCP, '_env', num2str(varEnv)])
  
    %% 6 - make plots
    get_plots(simu, i);
    
end