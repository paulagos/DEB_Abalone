%% Subfunction 
function [value,isterminal,direction] = event_spawning(t, EVHR, simu, fluxOut)

%% Unpack state vars
E    = EVHR(1);  % J, reserve energy
V    = EVHR(2);  % cm^3, structural volume
E_R  = EVHR(4);  % J, reproduction buffer

%% Load temperature conditions

if simu.envT == 1
    T = simu.Tinit;
else
    T = interp1(0:(length(simu.tempData)-1), simu.tempData, t, 'linear');
end

% Observable quantities
Ww_V      =  simu.par.d_V * V + (E + E_R) * simu.cPar.w_E / simu.par.mu_E;  % g, total wet-weight
G_Ww      =  E_R * simu.cPar.w_E / simu.par.mu_E;  % g, total wet-weight
GSI       =  G_Ww ./ (Ww_V - G_Ww);                          % -, Gonado somatic index (GSI)

value = ~(GSI >= 0.2 & T >= 17); % out the GSI and T threshold in the parameters

isterminal = [1]; % modified to stop
direction  = [0]; 
end