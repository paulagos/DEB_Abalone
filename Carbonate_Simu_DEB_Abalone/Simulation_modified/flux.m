function d = flux(t, EVHR, simu, par, cPar, fluxOut)
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

%% Unpack parameters
% Temperature correction parameters
T_ref = par.T_ref;
T_A   = par.T_A; 
T_L   = par.T_L;
T_AL  = par.T_AL;
T_H   = par.T_H;
T_AH  = par.T_AH;

% mobilisation, maintenance, growth & reproduction
kap   = par.kap;       %  -, allocation fraction to soma = growth + somatic maintenance
kap_R = par.kap_R;     %  -, reproduction efficiency
kap_X = par.kap_X;     %  -, digestion efficiency of food to reserve
p_M   = par.p_M;       % , J/d.cm^3, [p_M], vol-specific somatic maintenance
p_Am  = cPar.p_Am;
v     = par.v;

% life stages: E_H is the cumulated energy from reserve invested in maturation
E_Hb = par.E_Hb;      % J, E_H^b, maturity at birth
E_Hp = par.E_Hp;      % J, E_H^p, maturity at puberty
l_T  = cPar.l_T;
v_Hb = cPar.v_Hb;
v_Hj = cPar.v_Hj;
v_Hp = cPar.v_Hp;

% Additional parameters
K      = cPar.p_Am / (par.kap_X * par.F_m); % 22, same unit as food, half-saturation coefficient
g      = cPar.g;        % -, ﻿energy investment ratio
k_M    = cPar.k_M;      % 26
k_J    = par.k_J; 
E_G    = par.E_G;
d_V    = par.d_V;   
kap_LR = par.kap_LR;   % yield of reproduction tissue lysis
kap_LV = par.kap_LV;   % yield of structural  tissue lysis
rho_V  = par.rho_V ;   % J g^-1, energy content of 1g of structure 
k      = k_J / k_M;    % ratio of maturity and somatic maintenance rate coeff
L_m    = cPar.L_m;     % cm, maximun lenght
del_M  = par.del_M;    % -, shape coefficient for adults

% pCO2 correction parameters 
CM_H = par.CM_H;  % uatm, Higher pCO2 boundary for maintenace cost
CM_L = par.CM_L;  % uatm, Lower pCO2 boundary for maintenacne cost
q    = par.q;     % Jcm2 uatm-1 d-1, maintenace coefficient related to pCO2

% pars_tj = [g; k; l_T; v_Hb; v_Hj; v_Hp];
% [~, ~, ~, ~, l_p, ~, ~, ~, ~, ~] = get_tj(pars_tj, f);

% L_p  = L_m * l_p;                    % cm, structural length at puberty at f
% Lw_p = L_p / del_M;                 % cm, length at puberty at f
% V_p  = (del_M * Lw_p)^3;             % volume at puberty
V_p  = 2.1;             % volume at puberty
 
%% Load environmental conditions

if simu.envT == 1
    T = simu.Tinit;
else
    T = temp(t);
end
    
if simu.envX == 1
    X = simu.Xinit;
else
    X = food(t);
end

if simu.envpCO2 == 1
    pCO2 = simu.pCO2init;
else
    pCO2 = acid(t);
end

if simu.RCP == "RCP2"
    T = T + t * 0.64 / (365 * 10);
    pCO2 = pCO2 + t * 10 / (365 * 10);
end

if simu.RCP == "RCP3"
    T = T + t * 0.95 / (365 * 10);
    pCO2 = pCO2 + t * 15 / (365 * 10);
end

if simu.RCP == "RCP4"
    T = T + t * 1.5 / (365 * 10);
     pCO2 = pCO2 + t * 120 / (365 *10);
end

%% Unpack state vars
E    = EVHR(1);  % J, reserve energy
V    = EVHR(2);  % cm^3, structural volume
E_H  = EVHR(3);  % J, cumulated energy inversted into maturity
E_R  = EVHR(4);  % J, reproduction buffer

%% Functional response
f  = X / (X + K);

%% Temperature correction factors
s_A = exp(T_A/ T_ref - T_A ./ T);  % Arrhenius factor

s_L_ratio = (1 + exp(T_AL/ T_ref - T_AL/ T_L)) ./ ...
    (1 + exp(T_AL ./ T   - T_AL/ T_L));
s_H_ratio = (1 + exp(T_AH / T_H - T_AH/ T_ref)) ./ ...
    (1 + exp(T_AH / T_H - T_AH ./ T  ));

c_T = s_A .* ((T <= T_ref) .* s_L_ratio + (T > T_ref) .* s_H_ratio);

% Temperature corrected parameters
p_AmT = p_Am * c_T;                  % v, temp corrected energy conductance
p_XmT = p_AmT / par.kap_X;
k_JT  = k_J * c_T;

%% food quality modifier

kap_X = min(1, kap_X + simu.foodQ); 

%% DEB model Fluxes
if E_H < E_Hb
    pX = 0;                                                                        % J d-1, feeding
    pA = 0;                                                                        % J d-1, assimilation (eq. 2.3 p.35)
    pM = (p_M * V + q * min(max(pCO2 - CM_L, 0), max(CM_H - pCO2, 0)) * V^(2/3)) * c_T; % J d-1,  Somatic maintenance
    pJ = k_JT * E_H;                                                               % J d-1, Maturiry maintenance (eq.2.19 p.47) 
    pC = (E/V) * (E_G * v * V^(2/3) + p_M * V) / (kap * E/V + E_G) * c_T;          % J d-1, Reserve mobilization
    pG = kap * pC - pM;                                                            % J d-1, Growth
    pR = (1 - kap) * pC - pJ;                                                      % J d-1, Reproduction (eq.2.55 p.69)
    pD = pM + pJ + (1 - kap_R) * pR;                                               % J d-1, Dissipation
else
    pX = f * p_XmT * V^(2/3);                                                      % J d-1, feeding
    pA = kap_X * pX;                                                               % J d-1, Assimilation (eq.2.3 page 36)
    pM = (p_M * V + q * min(max(pCO2 - CM_L, 0), max(CM_H - pCO2, 0)) * V^(2/3)) * c_T; % J d-1, Somatic maintenance
    pJ = k_JT * E_H;                                                               % J d-1, Maturiry maintenance (eq.2.19 p.47) 
    pC = (E/V) * (E_G * v * V^(2/3) + p_M * V) / (kap * E/V + E_G) * c_T;          % J d-1, Reserve mobilization
    pG = kap * pC - pM;                                                            % J d-1, Growth
    pR = (1 - kap) * pC - pJ;                                                      % J d-1, Reproduction (eq.2.55 p.69)
    pD = pM + pJ + (1 - kap_R) * pR;                                               % J d-1,Dissipation
end

pLR    = 1 / kap_LR * (pM - kap * pC - ((1 - kap) * pC - pJ));
pLV    = 1 / kap_LV * (pM - kap * pC - ((1 - kap) * pC - pJ));
kap_pC = kap * pC;

%% Differential equations
 
dE = pA - pC; % dE/dt energy reserves
dV = pG / E_G; % dV/dt structural volume

if E_H < E_Hp
    dH = pR; % dEH/dt maturity
    dR = 0;  % dER/dt reproduction
    
else
    dH = 0;  % dEH/dt maturity
    dR = pR; % dR/dt reproduction
end

%
%     dE = pA - pC; % dE/dt energy reserves
%     
%     if kap * pC < pM % if in starvation, fraction from mobilization does not cover somatic maintenacce costs.
%  
%         dH = 0; % there is no maturation flux
%         
%         if V >= V_p && E_R >= pLR % if volume is larger that volume at puberty and reproduction buffer is larger than the amount of energy taken from the reproduction buffer to cover somatic maintenance
%             dV = 0;     % dV/dt structural volume
%             dR = -pLR;  % dR/dt reproduction buffer     
%         else
%             dV = -pLV / (rho_V * d_V);
%             dR = 0;
%         end
%         
%     else
%         dV = pG / E_G; % dV/dt structural volume
% 
%         if E_H < E_Hp
%             dH = pR; % dEH/dt maturity
%             dR = 0;  % dER/dt reproduction
%         else
%             dH = 0;  % dEH/dt maturity
%             dR = pR; % dR/dt reproduction
%         end
%     end
    
    if fluxOut
        d = [pX; pA; pM; pJ; pC; pG; pR; pD; pLR; pLV; kap_pC];
        
    else
        d = [dE; dV; dH; dR];
    end
