function d = flux(t, EVHR, simu, fluxOut)
%---------------------------------------------------------------
% Define differential equations of the state variables
%
%       t: n-vector with time points
% EVHR, 4-vector with state variables:
%         (1) E, J, reserve energy
%         (2) V, cm^3, structural volume
%         (3) E_H, J, cumulated energy inversted into maturity
%         (4) E_R, J, reproduction buffer
%         (5) E_S, cm^3, shell structural volume
%
% d: 4-vector with d/dt E, V, E_H, E_R
%
% called by : indiv.m,
% calls: food.m, temp.m,
%
% 2013/03/15 - Laure Pecquerie
% 2023/01/31 - Paulo F. Lagos
%--------------------------------------------------------------
 
%% Load environmental conditions
if simu.envT == 1
    T = simu.Tinit;
else
    T = interp1(0:(length(simu.tempData)-1), simu.tempData + 273.15, t, 'linear');
end
    
if simu.envX == 1
    X = simu.Xinit;
else
    X = food(t,simu);
end

if simu.envpCO2 == 1
    pCO2 = simu.pCO2init;
else
    pCO2 = interp1(0:(length(simu.pco2Data)-1), simu.pco2Data, t, 'linear');
end


%% Unpack state vars
E    = EVHR(1);  % J, reserve energy
V    = EVHR(2);  % cm^3, structural volume
E_H  = EVHR(3);  % J, cumulated energy inversted into maturity
E_R  = EVHR(4);  % J, reproduction buffer

%% Functional response
f  = X / (X + simu.K); 

%% aceleration factor
l_j = simu.cPar.v_Hj^(1/3);
l_b = simu.cPar.v_Hb^(1/3);
l = V.^(1/3) / simu.cPar.L_m;

if E_H < simu.par.E_Hb
    s_M = 1;
    
elseif simu.par.E_Hb <= E_H && E_H < simu.par.E_Hj
    s_M = l /l_b;
    
elseif E_H >= simu.par.E_Hj
    s_M = l_j / l_b;

else
    keyboard
end

%% Temperature correction factors
s_A = exp(simu.par.T_A / simu.par.T_ref - simu.par.T_A ./ T);  % Arrhenius factor

s_L_ratio = (1 + exp(simu.par.T_AL / simu.par.T_ref - simu.par.T_AL / simu.par.T_L)) ./ ...
    (1 + exp(simu.par.T_AL ./ T  - simu.par.T_AL / simu.par.T_L));
s_H_ratio = (1 + exp(simu.par.T_AH / simu.par.T_H - simu.par.T_AH/ simu.par.T_ref)) ./ ...
    (1 + exp(simu.par.T_AH / simu.par.T_H - simu.par.T_AH ./ T  ));

c_T = s_A .* ((T <= simu.par.T_ref) .* s_L_ratio + (T > simu.par.T_ref) .* s_H_ratio);

% Temperature corrected parameters
p_AmT = simu.cPar.p_Am * c_T;              % v, temp corrected energy conductance
k_JT  = simu.par.k_J * c_T;

%% food quality modifier
kap_X = min(1, simu.par.kap_X - simu.foodQ); 

%% applies zoom factors using primary scaling relationships
E_Hj  = simu.zoom.^3 * simu.par.E_Hj;
E_Hp  = simu.zoom.^3 * simu.par.E_Hp;
p_AmT = simu.zoom * p_AmT;

%% DEB model Fluxes
if E_H < E_Hj
    pX = f  * (s_M * p_AmT / simu.par.kap_X) * V^(2/3);                                    % feeding
    pA = kap_X * pX; % * c_T;                                                        % J d-1, assimilation (eq. 2.3 p.35)
    pM = (simu.par.p_M * V + simu.par.q_larva * min(max(pCO2 - simu.par.CML_larva, 0),...
        (simu.par.CMH_larva -  simu.par.CML_larva)) * V^(2/3)) * c_T; % somatic mnaintenace 
    pJ = k_JT * E_H;                                                              % J d-1, Maturiry maintenance (eq.2.19 p.47) 
    pC = (E / V) * ((simu.par.E_G * s_M * simu.par.v * V^(2/3) + simu.par.p_M * V) /...
        (simu.par.kap * E / V + simu.par.E_G)) * c_T; % Reserve mobilization
    pG = simu.par.kap * pC - pM;                                                           % J d-1, Growth
    pR = (1 - simu.par.kap) * pC - pJ;                                                     % J d-1, Reproduction (eq.2.55 p.69)
    pD = pM + pJ + (1 - simu.par.kap_R) * pR;                                              % J d-1, Dissipation
else
    pX = f * (s_M * p_AmT / simu.par.kap_X) * V^(2/3);                                                   % J d-1, feeding
    pA = kap_X * pX; % * c_T;                                                      % J d-1, Assimilation (eq.2.3 page 36)
    pM = (simu.par.p_M * V + simu.par.q_adult * min(max(pCO2 - simu.par.CML_adult, 0),...
        (simu.par.CMH_adult -  simu.par.CML_adult)) * V^(2/3)) * c_T;  % J d-1, Somatic maintenance
    pJ = k_JT * E_H;                                                            % J d-1, Maturiry maintenance (eq.2.19 p.47)
    pC = (E / V) * ((simu.par.E_G * s_M * simu.par.v * V^(2/3) + simu.par.p_M * V) /...
        (simu.par.kap * E / V + simu.par.E_G)) * c_T;       % J d-1, Reserve mobilization
    pG = simu.par.kap * pC - pM;                                                         % J d-1, Growth
    pR = (1 - simu.par.kap) * pC - pJ;                                                   % J d-1, Reproduction (eq.2.55 p.69)
    pD = pM + pJ + (1 - simu.par.kap_R) * pR;                                            % J d-1, Dissipation
end


%% Differential equations
 
dE = pA - pC; % dE/dt energy reserves

dV = (pG / simu.par.E_G); % dV/dt structural volume

if E_H < E_Hp %simu.par.E_Hp
    dH = (1 - simu.par.kap) * pC - k_JT * E_H;
    dR = 0;  % dER/dt reproduction
else
    dH = 0;   % dEH/dt maturity
    dR = simu.par.kap_R * (1 - simu.par.kap) * pC - k_JT * E_Hp; %simu.par.E_Hp;
end

if fluxOut
    d = [pX; pA; pM; pJ; pC; pG; pR; pD];
    
else
    d = [dE; dV; dH; dR];
    
end
   
