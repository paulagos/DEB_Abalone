%---------------------------------------------------------------
% Compute model predictions for the European abalone (Haliotis tuberculata)
%   with standard DEB equations (Kooijman 2010)
%
% using parameter estimates for the European abalone.
%
% From larva to adult stage under several climate scenarios and gergraphic locations along
% the coast of france.
% Model predictios include spawning event triggered by changes in sea suerface temperature (SST) temperature 
% and gonado somatic index (GSI)
% Options include predictions with intra-specific variation. (several individuals)
%
% simu : structure with individual features and model predictions
%
% calls :   init.m
%           indiv.m
%           get_fluxes.m
%           get_obs.m
%           get_plots.m
%           get_plots_zoom.m
%
% Additional funcions calls:
%
% 2013/03/15 - Laure Pecquerie
% 2023/03/12 - Paulo F. Lagos
%---------------------------------------------------------------

%%%%%%%%%%%%%  START OF OPTIONS TO BE SET BY USERS  %%%%%%%%%%%

%%%%%%%%%%%%%%     Change username accordingly      %%%%%%%%%%%%
%%%  and set path to folder containing results.mat file      %%%

username = 'paulo';

switch username
    
    case 'paulo'
        filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';
        results  = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/Simulation_DEB_Abalone/Simulation_files/DEB_model v.1.0/Results/';
    case 'laure'
        filename = 'to be defined by Laure'; 
        results  = 'to be defined by Laure'; 
    case 'jonathan'
        filename = 'to be defined by Jonathan';
        results  = 'to be defined by Jonathan';
    case 'sabine'
        filename = 'to be defined by Sabine';
        results  = 'to be defined by Sabine';
end

timeStamp = datestr(now,'yyyymmmdd');
saveDir = [results, 'simu_abal_', timeStamp, '/'];
if 7 ~= exist(saveDir, 'dir')
    mkdir(saveDir)
end


doPlotsSet1 = 1; % Set to 1 to plot SET 1,; set to 0 to plot SET 1 and SET  2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Set simulation Scenarios  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
SSP_scenarios = {'Current'}; %'Current','SSP1s19','SSP1s26','SSP1s45','SSP3sBaseline','SSP5sBaseline'};
all_locations = {'brest','roscoff','med'}; % locations: 'brest', 'roscoff', 'nantes','larochelle','med'
Food_quality  = [0]; % [0.3, 0.15, 0] percentages of decrease in food quality use any other combination;
variable_env  = [2]; % Constant (1) or Variable (2) environmental conditions
zoom_f        = [1]; % Number of idividuals to model) Range of zoom factors from 0.5 to 2 1 is default
gen_number    = [1]; % Increase number of simulations in times, For end of the century = 4
born_date     = [152]; % Year day whem animal begins life (1 is 1st of January, 365 is 31st December; any number in between is fine. we use  152 is 1st of June)
%
%%%%%%%%%% end of options for user %%%%%%%%%%%%%%%%%%%%%%%%%%
%
numberOfScenarios = length(SSP_scenarios)*length(all_locations)*...
    length(Food_quality)*length(variable_env)*length(zoom_f);
scenNo = 1;


close all

for h = 1:length(Food_quality)
    envFoodQ = Food_quality(h);
    
    for i = 1:length(SSP_scenarios)
        SSP = SSP_scenarios{i};
        
        for j = 1:length(variable_env)
            varEnv = variable_env(j);
            
            for k = 1:length(all_locations)
                location = all_locations{k};
                
                for m = 1:length(zoom_f)
                    zoom = zoom_f(m);
                    
                    sprintf('Running Scenario: %d of %d',scenNo,numberOfScenarios)
                    tic
                    
                    %% 1 - initialize time, parameters, etc
                    simu = init(filename, SSP, location, varEnv, envFoodQ, gen_number, zoom, born_date);%, pCO2, temp);
                    
                    for l = 1:gen_number
                        
                        if (gen_number > 1)
                            sprintf('Running Generation: %d of %d',l,gen_number)
                        end
                        
                        %% 2 - calculate state variables
                        simu.tEVHR(:,:,l) = indiv(filename, simu, l);
                        
                        %% 3 - calculate model fluxes
                        simu.tFluxes(:,:,l) = get_fluxes(filename, simu, l);
                        
                        %% 4 - calculate observable quantities
                        simu.obs(:,:,l) = get_obs(simu, l);
                        
                    end
                    
                    %% 5 - copy data structure to environmnet
                    simuSaveName = ['simu_foodQ', num2str(h), '_' ,SSP, '_',location ,'_env', num2str(varEnv), '_z', num2str(zoom*100)];
                    eval( [simuSaveName, '= simu;'] )
                    
                    %% 6 - save data structure to folder
                    resultFileName  = [saveDir, simuSaveName '.mat'];
                    save(resultFileName, 'simu')
                    
                    sprintf('Simulation finished for Scenario %d:',scenNo)
                    toc
                    scenNo = scenNo+1;
                end
            end
        end
    end
end

toc

if doPlotsSet1
    
    for h = 1:length(Food_quality)
        
        for i = 1:length(SSP_scenarios)
            SSP = SSP_scenarios{i};
            
            for j = 1:length(all_locations)
                location = all_locations{j};
                varEnv = 2; % for ploting  constant (1) or variable (2) environmental conditions
                
                % Dynamic evaluation to get the simu structures
                eval( ['simu = simu_foodQ', num2str(h), '_' ,SSP, '_',location,'_env', num2str(varEnv),'_z', num2str(zoom*100), ';'] )
                
                %% 6 - Export simulation outputs structure to a CSV file
                csvFileName = ['simu_foodQ', num2str(h), '_', SSP, '_', location, '_env', num2str(varEnv),'_z', num2str(zoom*100), '.csv'];
                exportStructToCSV(simu, csvFileName);
                
                %% 7 - make plots (set 1)
                get_plots(simu, h, i, gen_number);
                
                
            end
        end
    end
end

if any(zoom_f ~= 1)
    
    varEnv = 2; % for ploting  constant (1) or variable (2) environmental conditions
    for h = 1:length(Food_quality)
        
        for i = 1:length(SSP_scenarios)
            SSP = SSP_scenarios{i};
            
            for j = 1:length(all_locations)
                location = all_locations{j};
                
                for k = 1:length(zoom_f)
                    zoom = zoom_f(k);
                    
                    % Dynamic evaluation to get the simu structures
                    eval( ['simu = simu_foodQ', num2str(h), '_' ,SSP, '_',location,'_env', num2str(varEnv),'_z', num2str(zoom*100), ';'] )
                    
                    %% 7 - make plots for zoom factors (set 2)
                    get_plots_zoom(filename, simu, length(SSP_scenarios), length(all_locations), gen_number, j+(i-1)*length(all_locations), SSP, location);
                    
                end
                
            end
        end
    end
end



