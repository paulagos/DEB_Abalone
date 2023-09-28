function fluxes = get_fluxes(simu)

% unpack parameters values
filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/results_Haliotis_tuberculata2.mat';
res = load(filename);
par  = res.par;          % get primary parameters
cPar = parscomp_st(par); % get compound parameters

% get time vector
tVec = simu.tEVHRSO(:,1);

% set fluxOut to 1 to get fluxes rather than state variables from flux.m
fluxOut = 1;

% set up matrix to hold output
fluxes = zeros(length(tVec),11); % 10 is the number of fluxes output from flux.m
fluxes(:,1) = tVec;

for i = 1:length(tVec)
    
    t = tVec(i);
    EVHRSO = simu.tEVHRSO(i,2:size(simu.tEVHRSO,2));
    
    fluxes(i,2:size(fluxes,2)) = flux(t, EVHRSO, simu, par, cPar, fluxOut);
    
end