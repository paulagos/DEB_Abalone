function d = flux(t, EVHRSO, simu, par, cPar)
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

%% Load parameters

% mobilisation, maintenance, growth & reproduction
kap   = par.kap;       % 7, -, allocation fraction to soma = growth + somatic maintenance
kap_R = par.kap_R;     % 8, -, reproduction efficiency
p_M   = par.p_M;       % 9, J/d.cm^3, [p_M], vol-specific somatic maintenance
p_Am  = cPar.p_Am;

% life stages: E_H is the cumulated energy from reserve invested in maturation
E_Hb = par.E_Hb;      % 13, J, E_H^b, maturity at birth
E_Hp = par.E_Hp;      % 14, J, E_H^p, maturity at puberty

% Additional parameters
K      = cPar.p_Am / (par.kap_X * par.F_m); % 22, same unit as food, half-saturation coefficient
g      = cPar.g;        % 24
k_M    = cPar.k_M;      % 26
E_G    = par.E_G;
d_V    = par.d_V;
kap_LR = 0.79; % yield of reproduction tissue lysis
kap_LV = 0.79; % yield of structural  tissue lysis
rho_V   = 19600; % energy content of 1g of structure 

% pCO2 correction parameters (FOR DEVELPMENT ONLY, DELETE LATER!!!)
CX_L = 1010; % uatm,  Lower boundary for ingestion
CX_H = 6779; % uatm,  Upper boundary for ingestion
CM_H = 1300;  % uatm, Higher pCO2 boundary for maintenace cost
CM_L = 400;   % uatm, Lower pCO2 boundary for maintenacne cost
q    = 0.009; % Jcm2 uatm-1 d-1, maintenace coefficient related to pCO2

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

if simu.envCO2 == 1
    CO2 = simu.CO2init;
else
    CO2 = acid(t);
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

if CO2 <= CX_L
    c_X = 1;
end

if CX_L < CO2 && CO2 < CX_H
    c_X = (CX_H - CO2) / (CX_H - CX_L);
end

if CO2 >= CX_H
    c_X = 0;
end

% Temperature corrected parameters
p_AmT = c_T * p_Am;
v_T   = c_T * par.v;                   % v, temp corrected energy conductance
kT_M  = c_T * k_M;
p_XmT = p_AmT / par.kap_X;
p_MT  = c_T * p_M;
k_JT  = c_T * par.k_J;
p_TT  = c_T * par.p_T;  

% Test parameters (FOR DEVELPMENT ONLY, DELETE LATER!!!))
E_m   = p_AmT / par.v;              % Jcm3, maximun reserve density
e     = (E / V) / E_m;              % -,  scaled energy density
L_Vm  =  v_T / (g * kT_M);          % cm, maximum volumetric lenght
V_p   = 57; % For now this need to change later

%% DEB model Fluxes

if E_H < E_Hb
    pX = 0;                              % feeding
    pA = 0;                              % assimilation (eq. 2.3 p.35)
    pM = p_MT * V;                       % somatic mnaintenace
    pJ = k_JT * E_H;                     % Maturiry maintenance (eq.2.19 p.47)
    pT = p_TT * V^(2/3);                 % Surface-related somatic maintenance
    pS = pM + pT;                                                 % Somatic maintenance (eq.2.17 p.42)
    pC = (E/V) * (E_G * v_T * V^(2/3) + pS) / (kap * E/V + E_G ); %  Mobilisation (eq. 2.12 p.37)
    pG = kap * pC - pS;                                                                  % Growth
    pR = (1 - kap) * pC - pJ;                                      % Reproduction (eq.2.55 p.69)
    pD = 0;                                                        % dissipation
else
    pX = f * p_XmT * V^(2/3); 
    pA = par.kap_X * pX * c_X;                                                           % Assimilation (eq.2.3 page 36)
    pM = (p_M * V + q * min(max(CO2 - CM_L, 0), CM_H - CM_L * V^(2/3))) * c_T;            % Somatic maintenance
    pJ = min(V, V_p) * ((1 - kap) / kap) * p_M * c_T;                                    % Maturiry maintenance (eq.2.19 p.47)                                                      
    pC = (E/V) * ((par.E_G * par.v * V^(2/3) + p_M * V) / (kap * E/V + par.E_G)) * c_T;  % Reserve mobilization
    pG = kap * pC - pM;                                                                  % Growth
    pR = (1 - kap) * pC - pJ;                                                            % Reproduction (eq.2.55 p.69)
    pD = pM + pJ + (1 - kap_R) * pR;                                                     % Dissipation
end

pLR = 1/kap_LR * (pM - kap * pC - ((1 - kap) * pC - pJ)); 
pLV = 1/kap_LV * (pM - kap * pC - ((1 - kap) * pC - pJ));

% coefficients for pA, pD and pG (FOR DEVELPMENT ONLY, DELETE LATER!!!)
vA = 0.01;       % cm d-1, coupling coefficient to assimilation
vD = 2.37E-04;   % cm d-1, coupling coefficient to dissipation
vG = 0.1;        % cm d-1, coupling coefficient to growth
aA = 0;
aG = 1;
aD = 0;


%% Differential equations
%
    dE = pA - pC; % dE/dt energy reserves
    
    if kap * pC < pM % if in starvation condition fraction from mobilization does not cover somatic maintenacce costs.
        
        dS = (1 / cPar.p_Am) *  (vG * pG + vD * pD);              % dS/dt shell volume
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
        dV = pG / par.E_G; % dV/dt structural volume
        dS = (1 / cPar.p_Am) *  (vA * pA + vG * pG + vD * pD); % dS/dt shell volume
        dO = (aA * vA * pA + aG * vG * pG + aD * vD * pD) / (vA * pA + vG * pG + vD * pD); % dO/dt shell opacity (calcified structure)
        
        if E_H < E_Hp
            dH = pR; % dEH/dt maturity
            dR = 0;  % dER/dt reproduction
        else
            dH = 0;  % dEH/dt maturity
            dR = pR; % dR/dt reproduction
        end
    end

  d = [dE; dV; dH; dR; dS; dO];
