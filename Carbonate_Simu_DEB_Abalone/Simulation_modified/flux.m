function d = flux(t, EVHRSO, simu, par, cPar, fluxOut)
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
% l_T  = cPar.l_T;
% v_Hb = cPar.v_Hb;
% v_Hj = cPar.v_Hj;
% v_Hp = cPar.v_Hp;

% Additional parameters
K      = cPar.p_Am / (par.kap_X * par.F_m); % 22, same unit as food, half-saturation coefficient
g      = cPar.g;        % -, ﻿energy investment ratio
k_M    = cPar.k_M;      % 26
k_J    = par.k_J; 
E_G    = 2500;
d_V    = par.d_V;   
kap_LR = par.kap_LR;   % yield of reproduction tissue lysis
kap_LV = par.kap_LV;   % yield of structural  tissue lysis
rho_V  = par.rho_V ;   % J g^-1, energy content of 1g of structure 
k     = k_J / k_M;      % ratio of maturity and somatic maintenance rate coeff

% coefficients for pA, pD and pG (FOR DEVELPMENT ONLY, DELETE LATER!!!)
vA = 1.3;        % cm d-1, coupling coefficient to assimilation
vD = 1.2;        % cm d-1, coupling coefficient to dissipation
vG = 1.2;        % cm d-1, coupling coefficient to growth
aA = 0;
aG = 1;
aD = 0;

% pCO2 correction parameters (FOR DEVELPMENT ONLY, DELETE LATER!!!)
% CX_L = par.CX_L;  % uatm,  Lower boundary for ingestion
% CX_H = par.CX_H;  % uatm,  Upper boundary for ingestion
CM_H = par.CM_H;  % uatm, Higher pCO2 boundary for maintenace cost
CM_L = par.CM_L;  % uatm, Lower pCO2 boundary for maintenacne cost
q    = par.q;     % Jcm2 uatm-1 d-1, maintenace coefficient related to pCO2
 

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

E    = EVHRSO(1);  % J, reserve energy
V    = EVHRSO(2);  % cm^3, structural volume
E_H  = EVHRSO(3);  % J, cumulated energy inversted into maturity
E_R  = EVHRSO(4);  % J, reproduction buffer
E_S  = EVHRSO(5);  % Shell structural volume
O    = EVHRSO(6);  % shell opacity

%% Functional response

f  = X / (X + K);

%% Temperature and pCO2 correction factors

c_T   = exp(par.T_A / par.T_ref - par.T_A/ T);

% Temperature corrected parameters
p_AmT = p_Am * c_T;
v_T   = par.v * c_T;                   % v, temp corrected energy conductance
p_XmT = p_AmT / par.kap_X;
p_MT  = c_T * p_M;
k_JT  = par.k_J * c_T;
p_TT  = c_T * par.p_T;  

% Test parameters (FOR DEVELPMENT ONLY, DELETE LATER!!!))
V_p   = 1497; % For now this need to change later

%% DEB model Fluxes
if E_H < E_Hb
    pX = 0;                                                                        % feeding
    pA = 0;                                                                        % assimilation (eq. 2.3 p.35)
    pM = (p_M * V + q * min(max(pCO2 - CM_L, 0), max(CM_H - pCO2, 0)) * V^(2/3)) * c_T; % Somatic maintenance                                                            % Maturiry maintenance (eq.2.19 p.47)
    pJ = V * ((1 - kap) / kap) * p_M * c_T; 
    pC = (E/V) * (E_G * v * V^(2/3) + p_M * V) / (kap * E/V + E_G) * c_T;         % Mobilisation (eq. 2.12 p.37)
    pG = kap * pC - pM;                                                            % Growth
    pR = (1 - kap) * pC - pJ;                                                      % Reproduction (eq.2.55 p.69)
    pD = pM + pJ + (1 - kap_R) * pR;                                               % Dissipation
else
    pX = f * p_XmT * V^(2/3);                                                      % feeding
    pA = kap_X * pX;                                                               % Assimilation (eq.2.3 page 36)
    pM = (p_M * V + q * min(max(pCO2 - CM_L, 0), max(CM_H - pCO2, 0)) * V^(2/3)) * c_T; % Somatic maintenance
    pJ = V * ((1 - kap) / kap) * p_M * c_T;                              % Maturiry maintenance (eq.2.19 p.47) 
    pC = (E/V) * (E_G * v * V^(2/3) + p_M * V) / (kap * E/V + E_G) * c_T;        % Reserve mobilization
    pG = kap * pC - pM;                                                            % Growth
    pR = (1 - kap) * pC - pJ;                                                      % Reproduction (eq.2.55 p.69)
    pD = pM + pJ + (1 - kap_R) * pR;                                               % Dissipation
end

pLR = 1/kap_LR * (pM - kap * pC - ((1 - kap) * pC - pJ)); 
pLV = 1/kap_LV * (pM - kap * pC - ((1 - kap) * pC - pJ));

%% Differential equations
%
    dE = pA - pC; % dE/dt energy reserves
    
    if kap * pC < pM % if in starvation, fraction from mobilization does not cover somatic maintenacce costs.
        
        dS = vG * pG + vD * pD;                   % dS/dt shell volume
        dO = (aG *vG * pG + aD * vD * pD) / (vG * pG + vD * pD);  % dO/dt shell opacity (calcified structure)
        dH = 0; % there is no maturation flux
        
        if V >= V_p && E_R >= pLR % if volume is larger that volume at puberty and reproduction buffer is larger than the amount of energy taken from the reproduction buffer to cover somatic maintenance
            dV = 0;     % dV/dt structural volume
            dR = -pLR;  % dR/dt reproduction buffer     
        else
            dV = -pLV / (rho_V * d_V);
            dR = 0;
        end
        
    else
        dV = pG / E_G; % dV/dt structural volume
        dS = vA * pA + vG * pG + vD * pD; % dS/dt shell volume
        dO = (aA * vA * pA + aG * vG * pG + aD * vD * pD) / (vA * pA + vG * pG + vD * pD); % dO/dt shell opacity (calcified structure)

        
        if E_H < E_Hp
            dH = pR; % dEH/dt maturity
            dR = 0;  % dER/dt reproduction
        else
            dH = 0;  % dEH/dt maturity
            dR = pR; % dR/dt reproduction
        end
    end
    
    if fluxOut
  d = [pX; pA; pM; pJ; pC; pG; pR; pD; pLR; pLV];
    else
  d = [dE; dV; dH; dR; dS; dO];
    end
