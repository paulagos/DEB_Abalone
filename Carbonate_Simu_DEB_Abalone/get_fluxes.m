function fluxes = get_fluxes(filename, simu, gen)
%---------------------------------------------------------------
% unpack parameters values
res = load(filename);
par  = res.par;          % get primary parameters
cPar = parscomp_st(par); % get compound parameters
%---------------------------------------------------------------
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