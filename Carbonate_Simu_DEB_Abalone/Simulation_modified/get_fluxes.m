function fluxes = get_fluxes(simu, gen)

% % unpack parameters values
% filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';
% res = load(filename);
% par  = res.par;          % get primary parameters
% cPar = parscomp_st(par); % get compound parameters

% Prompt the user to select the .mat file
[filename, path] = uigetfile('*.mat', 'Select the .mat file');

if filename == 0
    % User canceled the selection
    disp('File selection canceled.');
    return;
end

% Construct the full file path
fullFilePath = fullfile(path, filename);

% Load the data
data = load(fullFilePath);

% Extract primary and compound parameters
par = data.par;
cPar = parscomp_st(par);

% Pull the variables into the workspace
vars_pull(par, cPar);
vars_pull(cPar);

% get time vector
tVec = simu.tEVHR(:, 1, gen);

% set fluxOut to 1 to get fluxes rather than state variables from flux.m
fluxOut = 1;

% set up matrix to hold output
fluxes = zeros(length(tVec),12); % 10 is the number of fluxes output from flux.m
fluxes(:,1) = tVec;

for i = 1:length(tVec)
    
    t = tVec(i);
    EVHR = simu.tEVHR(i,2:size(simu.tEVHR,2), gen);
    
    fluxes(i,2:size(fluxes,2)) = flux(t, EVHR, simu, par, cPar, fluxOut);
    
end