function [prdData, info] = predict_Haliotis_tuberculata(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
  % Filters for allowable parameters of the standard DEB model (std)
  filterChecks = k * v_Hp >= f^3 || ...         % constraint required for reaching puberty with f
                 ~reach_birth(g, k, v_Hb, f);   
             
  if filterChecks  
    info = 0;
    prdData = {};
    return;
  end  
  
  % compute temperature correction factors
  pars_T     = [T_A; T_L; T_H; T_AL; T_AH];
  TC_ab      = tempcorr(temp.ab, T_ref, pars_T);
  TC_tp      = tempcorr(temp.tp, T_ref, pars_T);
  TC_am      = tempcorr(temp.am, T_ref, pars_T);
  TC_tj      = tempcorr(temp.tj, T_ref, pars_T);
  TC_tW      = tempcorr(temp.tW, T_ref, pars_T);
  TC_LN      = tempcorr(temp.LN, T_ref, pars_T);
  TC_Ri      = tempcorr(temp.Ri, T_ref, pars_T);
  TC_avignon = tempcorr(temp.tW_Avignon, T_ref, pars_T);
  TC_wessel  = tempcorr(temp.tL_wessel, T_ref, pars_T);
  TC_roussel = tempcorr(temp.LGw_ee, T_ref, pars_T);
  TC_LW0     = tempcorr(temp.Lw0, T_ref, pars_T);
  
  %  parameter sets
  pars_tj = [g; k; l_T; v_Hb; v_Hj; v_Hp];
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  
  pars_R = [kap, kap_R, g, k_J, k_M, L_T, v, U_Hb, U_Hj, U_Hp];
  pars_UE0 = [V_Hb; g; k_J; k_M; v];

  %% Zero-variate data

  % initial
  U_E0  = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
  p_TAm = p_Am * TC_LW0;                       %, J/d.cm^2, {p_Am}, spec. temperature correted assimilation flux
  E_0   = U_E0 * p_TAm;                        % J, initial energy 
  Wd_0  = E_0 * w_E/ mu_E;                     % g, egg dry weight
  V_egg = E_0 * w_E / (mu_E * d_E) / del_Y^3;  % cm^3, volume of egg, From: Augustine, Litvak & Kooijman,(2011)
  Lw_0  = V_egg^(1/3);                         % cm, diameter of egg
  
  % life cycle
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  if info ~= 1 % numerical procedure failed
    fprintf('warning: invalid parameter value combination for get_tj \n')
    prdData = {};
    return
  end
  
  % birth
  aT_b = t_b / k_M / TC_ab;           % d, age at birth at f and T
  L_b = L_m * l_b;                    % cm, structural length at birth at f (post-metamorphosis)
  Lw_b = L_b / del_L;                 % cm, total length at birth at f (post-metamorphosis)
  V_b = (L_b / del_M).^3;             % cm^3, volume at birth (no pH) (post-metamorphosis)
  E_b = f * E_m * L_b^3;              % J, energy reserves at birth (no pH) (post-metamorphosis)
%   Ww_b = L_b^3 * (1 + f * w);        % g, wet-weight at birth
%   Wd_b = L_b^3 * d_V * (1 + f * w);  % g, dry-weight at birth
  
  % metamorphosis
  tT_j = t_j / k_M / TC_tj;           % d, age at metamorphosis
  L_j = L_m * l_j;                    % cm, structural length at metam
  Lw_j = L_j / del_L;                 % cm, total lenght at metamorphosis at f
%   Ww_j = L_j^3 * (1 + w * g);        % g, wet weight at metam
%   Wd_j = L_j^3 * d_V * (1 + f * w);  % g, dry weight at metamorphosis
    
  % puberty 
  tT_p = t_p  / k_M / TC_tp;          % d, age at puberty at f and T
  L_p = L_m * l_p;                    % cm, structural length at puberty at f
  Lw_p = L_p / del_M;                 % cm, length at puberty at f
  V_p = (del_M * Lw_p)^3;             % volume at puberty 
%   Wd_p = L_p^3 * d_V * (1 + f * w);  % g, dry weight at puberty 
  
   % ultimate
  L_i =  L_m - l_i;                    % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                   % cm, ultimate total length at f
  Ww_i = L_i^3 * (1 + f * w);          % g, ultimate wet-weight
 % Wd_i = L_i^3 * d_V * (1 + f * w);    % g, ultimate dry weight 
 % Ww_i = L_i^3 * (1 + f * w);          % g, shell-free wet weight
   
  % life span
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
  % reproduction
  RT_i = 365 * TC_Ri * reprod_rate_j(L_i, f, pars_R);   % #/d, ultimate reproduction rate at T
    
  % pack to output
  prdData.ab  = aT_b;
  prdData.tj  = tT_j;
  prdData.tp  = tT_p;
  prdData.am  = aT_m;
  prdData.Lw0 = Lw_0;
  prdData.Lb  = Lw_b;
  prdData.Lj  = Lw_j;
  prdData.Lp  = Lw_p;
  prdData.Li  = Lw_i;
%   prdData.Wdb = Wd_b;
%   prdData.Wdj = Wd_j;
%   prdData.Wdp = Wd_p;
%   prdData.Wwb = Ww_b;
%   prdData.Wwj = Ww_j;
  prdData.Wm = Ww_i;
%   prdData.Wdi = Wd_i;
  prdData.Ri  = RT_i;
  
  %% Uni-variate data
  % Basu (1997) time (d) ~ weight(g) (adults)
  kT_M  = k_M * TC_tW;
  rT_B  = rho_B * kT_M;
  tT_j  = (t_j - t_b)/ kT_M;                                    % time at metamorphosis (end of acceleration)
  L     = L_i - (L_i - L_j) * exp(-rT_B * (tW(:,1) - tT_j));    % cm, strucutural length
  EWw   = L.^3 * (1 + f * w);                                   % g,  total wet-weight
  EWw   = EWw + 0.5 * 365 * TC_tW * reprod_rate_j(L, f, pars_R) * Wd_0 / d_E; % g, total wet-weight including reproduction buffer
  
  % Bilbao (2010) Lenght (cm) ~ wet-weight (g)
  EWw_f = (LWw_f(:,1) * del_M).^3 * (1 + f * w);                 % g, expected total wet-weights females (shell included)
  EWw_m = (LWw_m(:,1) * del_M).^3 * (1 + f * w);                 % g, expected total wet-weight males (shell included)
  
  % Roussel (2019) Length (cm) ~ wet-weight (g)
  ELWw_roussel = (LWw_roussel(:,1) * del_M).^3 * (1 + f * w);     % g, expected total wet-weights (shell included)
  
  % Bilbao (2010) Total length (cm) ~ Batch fecundity (Females)
  EN = 365 * TC_LN * reprod_rate_j(LN(:,1), f, pars_R);          % # eggs, batch fecundity
  
  %!!!! PROBLEMS WITH STARTING CONDTIONS!!! aprox. 600 days old adults abalone, 
  %!!!! but, I am using starting conditions for the ODE for birth (not sure how to find other starting conditions that match the abalone age or extract the growth data from the ODE for a time interval that match the data)
  % Avignon (2020)  time(d) ~ wet-weight (g)
  EVHRS_1 = [E_b; V_b; E_Hb; 0; 0]; % initial conditions for state variables E,V,H,R and S
  [t, EVHRS_pH] = ode45(@dget_EVHRS_pH, [0 ; data.tL_Avignon(:,1)], EVHRS_1, [], p_Am, v, V_b, V_p, d_V, kap, kap_X, kap_R, p_M, k_J, p_T, E_G, E_Hb, E_Hp, q, CM_L, CM_H, kap_LR, kap_LV, rho_V, vA, vG, vD, T_ref, T_A, temp.tW_Avignon, pCO2.tW_Avignon, f); %[] in third no ODE options
  EVHRS_pH(1,:) = []; E = EVHRS_pH(:,1); V = EVHRS_pH(:,2); R = EVHRS_pH(:,4);
  EtW_Avignon = d_V .* V + (E + R) .* w_E / mu_E;                % g,  Wet-weight (including gonad)
  
  %!!!! PROBLEMS WITH STARTING CONDTIONS!!! aprox. 600 days old adults abalone, but I am using starting conditions for the ODE for birth
  % Avignon (2020)  time(d) ~ length (cm)
  EVHRS_1 = [E_b; V_b; E_Hb; 0; 0]; % initial conditions for state variables E,V,H,R and S
  [t, EVHRS_pH] = ode45(@dget_EVHRS_pH, [0 ; data.tL_Avignon(:,1)], EVHRS_1, [], p_Am, v, V_b, V_p, d_V, kap, kap_X, kap_R, p_M, k_J, p_T, E_G, E_Hb, E_Hp, q, CM_L, CM_H, kap_LR, kap_LV, rho_V, vA, vG, vD, T_ref, T_A, temp.tL_Avignon, pCO2.tL_Avignon, f); %[] in third no ODE options
  EVHRS_pH(1,:) = []; E = EVHRS_pH(:,1); V = EVHRS_pH(:,2);
  EtL_Avignon = V.^(1/3) ./ del_M;                                % cm, total length
  
  % Roussel (2019) lenght(cm) ~ gonad weight(g) ~ wet-weight (g) at f1
  ELGW_sl   = TC_roussel  * 365 * reprod_rate_j(LGw_sl(:,1), f_sl, pars_R) * Wd_0 / d_E; % g, Expected total gonad wet-weight at diet of Saccorhiza polyschides
  EWw_sl = (LWw_sl(:,1) * del_M).^3 * (1 + f_sl * w);            % g, expected total wet-weight at diet of Saccorhiza polyschides
  
  % Roussel (2019) lenght(cm) ~ gonad weight(g) ~ wet-weight (g) at f2
  ELGw_pp   = TC_roussel  * 365 * reprod_rate_j(LGw_pp(:,1), f_pp, pars_R) * Wd_0 / d_E; % g, Expected total gonad wet-weight at diet of Palmaria palmata
  EWw_pp = (LWw_pp(:,1) * del_M).^3 * (1 + f_pp * w);            % g, expected total wet-weight at diet of Palmaria palmata
  
  % Roussel (2019) lenght(cm) ~ gonad weight(g) ~ wet-weight (g) at f3
  ELGw_ee = TC_roussel  * 365 * reprod_rate_j(LGw_ee(:,1), f_ee, pars_R) * Wd_0 / d_E; % g, Expected total gonad wet-weight at diet of Enteromoepha sp.
  EWw_ee = (LWw_ee(:,1) * del_M).^3 * (1 + f_pp * w);            % g, expected total wet-weight at diet of Enteromoepha sp.
  
  % !!!!This one I can solve it!! still working on it thought, I just need to match the temperature interpolation output to the lenght of the data
  % Auzoux-Bordenave (2021)  time (days) ~ total length (cm)  
%   EVHRS_1 = [E_b; V_b; E_Hb; 0; 0]; % initial conditions for state variables E,V,H,R and S
%   [t, EVHRS_pH] = ode45(@dget_EVHRS_pH, [0 ; data.tL_Auzoux(:,1)], EVHRS_1, [], p_Am, v, V_b, V_p, d_V, kap, kap_X, kap_R, p_M, k_J, p_T, E_G, E_Hb, E_Hp, q, CM_L, CM_H, kap_LR, kap_LV, rho_V, vA, vG, vD, T_ref, T_A, temp.tL_Auzoux, pCO2.tL_Auzoux, f); %[] in third no ODE options
%   EVHRS_pH(1,:) = []; E = EVHRS_pH(:,1); V = EVHRS_pH(:,2);
%   EtL_Auzoux = V^(1/3) ./ del_M;                                    % cm, total length
  
  % Wessel (2018) time (dpf) - larval lenght (cm)
  kT_M  = k_M * TC_wessel;
  rT_j  = rho_j * kT_M;                                             % 1/d, growth rate larva
  L_bj = L_b * exp((tL_wessel(:,1) - aT_b) * rT_j/ 3);
  %L     = L_i - (L_i - L_j) * exp(-rT_B * (t(t >= tT_j) - tT_j));   % cm, expected length at time % cm, strucutural length
  EtL_wessel = L_bj ./ del_L;                                    % g, expected wet-weight
  
  % Chauvaud (2019) length (cm) ~ total weight (g)
  EWw_Chauvaud = (LW_chauvaud(:,1) * del_M).^3 * (1 + f * w);      % g, expected total wet-weights (shell included)
  
  % pack to output
  prdData.tW    = EWw;
  prdData.LWw_f = EWw_f;
  prdData.LWw_m = EWw_m;
  prdData.LWw_roussel = ELWw_roussel;
  prdData.LN    = EN;
  prdData.tL_Avignon  = EtL_Avignon;
  prdData.tW_Avignon  = EtW_Avignon;
  prdData.LGw_sl = ELGW_sl;
  prdData.LWw_sl = EWw_sl;
  prdData.LGw_pp = ELGw_pp;
  prdData.LWw_pp = EWw_pp;
  prdData.LGw_ee = ELGw_ee;
  prdData.LWw_ee = EWw_ee;
%   prdData.tL_Auzoux = EtL_Auzoux;
  prdData.tL_wessel = EtL_wessel;
  prdData.LW_chauvaud = EWw_Chauvaud;
  
end

%% Subfunction
function dEVHRS = dget_EVHRS_pH(~, EVHRS, p_Am, v, V_b,V_p, d_V, kap, kap_X, kap_R, p_M, ~, ~, E_G, E_Hb, E_Hp, q, CM_L, CM_H, kap_LR, kap_LV, rho_V , vA, vG, vD, T_ref, T_A, tT, pCO2, f)

  %  change in state variables during juvenile and adult stage at given temperatues and pH levels
  %  dELHR = dget_ELHR_pj(t, ELH)
  %  ELHR: 4-vector
  %   E: reserve E in (J)
  %   V: structural volume in(cm3)
  %   E_H: maturity E_H in (J)
  %   E_R: reproduction buffer E_R in (J)
  %   S: shell volume (cm3)
  %  dELHR: change in reserve, length, maturity and reproduction buffer
  
  %  unpack variables
  E = EVHRS(1); V = EVHRS(2); E_H = EVHRS(3); E_R = EVHRS(4); S = EVHRS(5);
  
  % Interpolation for temperatures
  if length(tT) > 1
      TC = tempcorr(interp1(tT(:,1),tT(:,2), t), T_ref, T_A);
  else
      TC = tempcorr(tT, T_ref, T_A); % -, temperature correction factor
  end
  
% Temperature corrected parameters
p_AmT = TC * p_Am;                
p_XmT = p_AmT / kap_X;

% DEB Fluxes
if E_H < E_Hb % egg and Larvae
    pX = 0;                                                                           % feeding
    pA = 0;                                                                           % assimilation (eq. 2.3 p.35)
    pM = (p_M * V + q * min(max(pCO2 - CM_L, 0), max(CM_H - pCO2, 0))) * V^(2/3) * TC; % somatic mnaintenace
    pJ = min(V, V_b) * ((1 - kap) / kap) * p_M * TC;                               % Maturiry maintenance (eq.2.19 p.47)
    pC = (E/V) * (E_G * v * V^(2/3) + p_M * V) / (kap * E/V + E_G ) * TC;          % Mobilisation (eq. 2.12 p.37)
    pG = kap * pC - pM;                                                               % Growth
    pR = (1 - kap) * pC - pJ;                                                         % Reproduction (eq.2.55 p.69)
    pD = pM + pJ + (1 - kap_R) * pR;                                                  % Dissipation
    
else % Post-metamorphosis, Juvenile and adults
    pX = f * p_XmT * V^(2/3); 
    pA = kap_X * pX * TC;                                                           % Assimilation (eq.2.3 page 36)
    pM = (p_M * V + q * min(max(pCO2 - CM_L, 0), max(CM_H - pCO2, 0))) * V^(2/3) * TC; % somatic mnaintenace
    pJ = min(V, V_p) * ((1 - kap) / kap) * p_M * TC;                                % Maturiry maintenance (eq.2.19 p.47)                                                      
    pC = (E/V) * ((E_G * v * V^(2/3) + p_M * V) / (kap * E/V + E_G)) * TC;             % Reserve mobilization
    pG = kap * pC - pM;                                                              % Growth
    pR = (1 - kap) * pC - pJ;                                                        % Reproduction (eq.2.55 p.69)
    pD = pM + pJ + (1 - kap_R) * pR;                                                 % Dissipation
end

pLR = 1/ kap_LR * (pM - kap * pC - ((1 - kap) * pC - pJ)); 
pLV = 1/ kap_LV * (pM - kap * pC - ((1 - kap) * pC - pJ));

% Differential equations

    dE = pA - pC; % dE/dt energy reserves
    
    if kap * pC < pM % if in starvation conditions, fraction from mobilization does not cover somatic maintenacce costs.
        
        dS = (1 / p_Am) *  (vG * pG + vD * pD);              % dS/dt shell volume
        dH = 0; % there is no maturation flux
        
        if V >= V_p && E_R >= pLR % if volume is larger that volume at puberty and reproduction buffer is larger than the amount of energy taken from the reproduction buffer to cover somatic maintenance
            dV = 0;     % dV/dt structural volume (or there is no growth)
            dR = -pLR;  % dR/dt reproduction buffer     
        else
            dV = -pLV / (rho_V * d_V);
            dR = 0;
        end
        
    else
        dV = pG / E_G; % dV/dt structural volume
        dS = (1 / p_Am) *  (vA * pA + vG * pG + vD * pD); % dS/dt shell volume
        
        if E_H < E_Hp
            dH = pR; % dEH/dt maturity
            dR = 0;  % dER/dt reproduction
        else
            dH = 0;  % dEH/dt maturity
            dR = pR; % dR/dt reproduction
        end
    end

  % pack derivatives
  dEVHRS = [dE; dV; dH; dR; dS];
end
