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
%           food.m 
%           temp.m
%           acid.m
%           get_obs.m 
%           get_plots 
%
% 2013/03/15 - Laure Pecquerie
% 2023/03/12 - Paulo F. Lagos
%---------------------------------------------------------------  

  clear all
  %close all
  
  %% 1 - initialize time, parameters, etc
  simu = init;
  
  %% 2 - calculate state variables
  simu.tEVHRSO = indiv(simu);
  
  %% 3 - calculate observable quantities
  simu.obs = get_obs(simu);
  
  %% 4 - make plots
    get_plots(simu);
  
  