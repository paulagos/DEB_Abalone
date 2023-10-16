function tEVHR = indiv(simu, gen)
%---------------------------------------------------------------
% Compute model predictions (numerical integration)
%
% simu: 1-structure with individual features (parameters, env...),;
%             see init.m
%
% tEVHRS : (nt,4) matrix with time and state variables
%         t, d,      time
%         E , J,     reserve energy
%         V , cm^3,  structural volume
%         E_H , J ,  cumulated energy inversted into maturity
%         E_R , J,   reproduction buffer
%         E_S, cm^3, shell structural volume
%         O  ,-,     shell opacity
%currently
% called by: main.m
% calls    :flux.m
%
% 2013/03/15 - Laure Pecquerie
% 2023/02/22 - Paulo F. Lagos
%--------------------------------------------------------------

%% unpack parameters values
filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';
res = load(filename);
par  = res.par;          % get primary parameters
cPar = parscomp_st(par); % get compound parameters
vars_pull(par);
vars_pull(cPar);

%%
tc = simu.t0 + (gen - 1) * simu.tm; % set current time to initial time
EVHR_tc = simu.EVHR_init; % set initial values of state variables
tEVHR = zeros(simu.tm+1, 5); % initialize outputs
gen_tm = gen * simu.tm;    % set the maximun simulation time span for each generation number
gen_ts = simu.ts + tc;

tic
i = 0; % initialize year index

while tc < gen_tm
    
    i = i + 0.5; %spawning twice every year
    tnext = gen_ts + 365 * i; % integrate between two spawning dates (two spawning each year)
    
    if tnext >= gen_tm
        tnext = gen_tm;
    end
    
    t = unique([tc, ceil(tc): floor(tnext)])';
    
    fluxOut = 0; % make flux output the usual state variable output, and not the fluxes
    [t, EVHR] = ode45(@(t, EVHR) flux(t, EVHR, simu, par, cPar, fluxOut), t, EVHR_tc); % integrate
    
    % indexing time for next generation
    tEVHR( (t - (gen - 1) * simu.tm +1) ,:) = [t, EVHR];
    
    % setting mortality
    %     V_diff = tEVHR((t(end)-(gen-1) * simu.tm - 180+1):(t(end)-(gen-1) * simu.tm),3) - tEVHR((t(end)-(gen-1) * simu.tm - 180):(t(end)-(gen-1)*simu.tm - 1),3);
    
    E_Hc = EVHR(end, 3);
    
    if E_Hc >= E_Hp %; E_Hp
        EVHR(end,4) = 0;
    end
    EVHR_tc = EVHR(end,:)';
    tc = t(end)+1;
    
    %     if sum(V_diff < 0) == 180 % if the animal has been starving for every day of the last 180 days, then it should die
    %         EVHR(end,4) = 0;
    
    %     end
    %
    
end
toc
end