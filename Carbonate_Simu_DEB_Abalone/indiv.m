function tEVHR = indiv(filename, simu, gen)
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
res = load(filename);
par  = res.par;          % get primary parameters
cPar = parscomp_st(par); % get compound parameters
vars_pull(par);
vars_pull(cPar);

%%
tc = simu.t0 + (gen - 1) * simu.tm; % set current time to initial time
EVHR_tc = simu.EVHR_init; % set initial values of state variables
tEVHR   = zeros(simu.tm+1, 5); % initialize outputs
gen_tm  = gen * simu.tm;    % set the maximun simulation time span for each generation number
tSpan   = [tc:gen_tm];

while tc < gen_tm
    

    %
    fluxOut = 0; % make flux output the usual state variable output, and not the fluxes
    
    options = odeset('Events', @event_spawning, 'MaxStep', 1);
 
    [tOut, EVHR, te,ye,ie] = ode45(@flux, tSpan, EVHR_tc, options, simu, fluxOut); % integrate
    
    if ( ~isempty(te) )
        
        % Check if we ended on non-integer and, if so, trim it out
        if ( round(te) ~= te )
            tOut = tOut(1:(end-1));
            EVHR = EVHR(1:(end-1),:);
        end
        
    end
    
    % saving output in predefined matrix
    tEVHR( (tOut - (gen - 1) * simu.tm +1) ,:) = [tOut, EVHR];
    
    if ( tOut(end) >= gen_tm )
        tc = tOut(end);
    elseif ( tOut(end) >= gen_tm-2 )
        
        nrows = gen_tm-tOut(end);
        tEVHR(2+simu.tm-nrows:1+simu.tm,2:end) = repmat(EVHR(end,:),nrows,1);
        tEVHR(2+simu.tm-nrows:1+simu.tm,1) = 1+gen_tm-nrows:gen_tm;
        tEVHR( (tOut(end)-(gen-1)*simu.tm+2), 4) = 0;
        if nrows == 2
            tEVHR((tOut(end)-(gen-1)*simu.tm+3), 4) = 0.001;
        end
        tc = gen_tm;
    else
        % Define new tSpan and initial EVHR for the sim
        tSpan = (tOut(end)+1):gen_tm;
        
        EVHR_tc = EVHR(end,:)';
        
        
        % Implement spawning, i.e. reset E_R to 0
        E_Hc = EVHR_tc(3);
        E_Hp = simu.zoom.^3 * simu.par.E_Hp;
        propRemaining = 0.2; % proportion of gonad remaining after spawning
        if E_Hc >= E_Hp
            EVHR_tc(4) = propRemaining * EVHR_tc(4);
        end
        
        tc      = tSpan(1);
    end
    
    
end

end