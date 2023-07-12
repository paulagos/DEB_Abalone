close all; 
global pets 

pets = {'Haliotis_tuberculata'}; 
check_my_pet(pets); 

estim_options('default'); 
estim_options('max_step_number',1000); 
estim_options('max_fun_evals',4000); 

%sb - symmetric bounded (default)
%su - symmetric unbounded
estim_options('loss_function', 'sb')

% pars_init_method':
%   0 - get initial estimates from automatized computation (default)
%   1 - read initial estimates from .mat file (for continuation)
%   2 - read initial estimates from pars_init file
estim_options('pars_init_method', 1);

%  'results_output':
%   0     - only saves data to .mat (no printing to html or screen and no figures) - use this for (automatic) continuations
%   1, -1 - no saving to .mat file, prints results to html (1) or screen (-1), shows figures but does not save them
%   2, -2 - saves to .mat file, prints results to html (2) or screen (-2), shows figures but does not save them
%   3, -3 - like 2 (or -2), but also prints graphs to .png files (default is 3)
%   4, -4 - like 3 (or -3), but also prints html with implied traits
%   5, -5 - like 4 (or -4), but includes related species in the implied traits
%   6     - like 5, but also prints html with population traits
estim_options('results_output', 2);

% 'method':
%   'nm' - use Nelder-Mead method;
%   'no' - do not estimate; 
estim_options('method', 'nm'); 

% for implied properties of the model
% load('results_Nyctiphanes_australis.mat');
% prt_report_my_pet({par, metaPar, txtPar, metaData}),

tStart = tic;
estim_pars; 
%main;       % Run the simulation code 
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));


% Save results to working folder
figures = dir('*.png');
html    = dir('*.html');
results = dir('*.mat');

results_file = {figures.name,...
    html.name,...
    results.name,...
    };
timeStamp = char(datetime('now'));
saveDir   = ['results_', timeStamp, '/'];
subDir    = [saveDir,'Trial 1_'];
mkdir(subDir);

for i = 1:length(results_file)
    copyfile(results_file{i},[subDir, results_file{i}])
end

