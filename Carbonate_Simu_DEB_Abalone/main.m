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

username = 'paulo';
%username = 'laure';

switch username
    case 'paulo'
        filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';
    case 'laure'
        filename = 'to be defined by Laure';
    case 'jonathan'
        filename = 'to be defined by Jonathan';
end

close all

SSP_scenarios = {'Current','SSP1_19','SSP1_26','SSP1_45','SSP3_baseline','SSP5_baseline'}; 
all_locations = {'marel'}; % add the other locations once the data has been cleaned 
Food_quality  = [0.3, 0]; % [0.3, 0.15, 0];
variable_env  = (2);  % constant (1) or variable (2) environmental conditions
gen_number    = [4];

for h = 1:length(Food_quality)
    envFoodQ = Food_quality(h);
    
    for i = 1:length(SSP_scenarios)
        SSP = SSP_scenarios{i};
        
        for j = 1:length(variable_env)
            varEnv = variable_env(j);
            
            for k = 1:length(all_locations)
                location = all_locations{k};
            
            %% 1 - initialize time, parameters, etc
            simu = init(filename, SSP, location, varEnv, envFoodQ, gen_number);%, pCO2, temp);
            
            for l = 1:gen_number
                
                %% 2 - calculate state variables
                simu.tEVHR(:,:,l) = indiv(filename, simu, l);
                
                %% 3 - calculate model fluxes
                simu.tFluxes(:,:,l) = get_fluxes(filename, simu, l);
                
                %% 4 - calculate observable quantities
                simu.obs(:,:,l) = get_obs(filename, simu, l);
                
            end
            
            %% 5 - save data structure
            eval( ['simu_foodQ', num2str(h), '_' ,SSP, '_',location ,'_env', num2str(varEnv), '= simu;'] )
        
            end
        end
    end
end

for h = 1:length(Food_quality)
    
    for i = 1:length(SSP_scenarios)
        SSP = SSP_scenarios{i};
        
        for j = 1:length(all_locations)
            location = all_locations{j};
            
            varEnv = 2; % for ploting  constant (1) or variable (2) environmental conditions
            
            eval( ['simu = simu_foodQ', num2str(h), '_' ,SSP, '_',location,'_env', num2str(varEnv), ';'] )
            
            %% 6 - make plots
            get_plots(filename, simu, i, h, gen_number);
          
           
        end
    end
end



