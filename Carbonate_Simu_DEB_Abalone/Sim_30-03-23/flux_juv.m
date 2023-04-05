function d = flux_juv(t, EVHRSO, simu)
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
  % Environmental conditions
  if simu.envT == 1
      T = simu.Tinit;   
  else
      T = temp(t); 
  end

  if simu.envX == 1
      X = simu.Xinit ;
  else
       X = food(t);
  end
  
   if simu.envCO2 == 1
      CO2 = simu.CO2init;
  else
      CO2 = acid(t);
  end
   
  %% unpack state vars
  E    = EVHRSO(1);  % J, reserve energy
  V    = EVHRSO(2);  % cm^3, structural volume
  E_H  = EVHRSO(3);  % J , cumulated energy inversted into maturity 
  E_R  = EVHRSO(4);  % J, reproduction buffer
  E_S  = EVHRSO(5);  % Shell structural volume
  O    = EVHRSO(6);  % shell opacity
  
  %% Read parameter values
  par   = simu.par;
  T1    = par(1);   % K, Reference temperature  
  TA    = par(2);   % K, Arrhenius temperature 
  kap_R = par(3);   % -, reproduction efficiency
  kap_X = par(4);   % -, digestion efficiency of food to reserve
  p_Am  = par(5);   % J/cm^2/d, {p_Am}, maximum surface-specific assimilation rate 
  v     = par(6);   % cm/d, energy conductance
  kap   = par(7);   % -, allocation fraction to soma = growth + somatic maintenance
  p_M   = par(9);   % J/d.cm^3, [p_M], vol-specific somatic maintenance
  p_T   = par(10);  % J/d.cm^2, {p_T}, surface-specific som maintenance
  k_J   = par(11);  % 1/d, maturity maint rate coefficient
  E_G   = par(12);  % J/cm^3, [E_G], spec cost for structure
  E_Hb  = par(13);  % J, E_H^b, maturity at birth
  E_Hp  = par(14);  % J, E_H^p, maturity at puberty
  K     = par(22);  % same unit as food, half saturation constant
  g     = par(24);  % -, energy investment ratio
  k_M   = par(26);  % 
  
  % scaled functional response
  f = X / (X + K); % -,   scaled functional response\
  
  % simplest temperature correction function
  c_T   = exp(TA/ T1 - TA/ T);
  
  % pCO2 correction function (FOR DEVELPMENT ONLY, DELETE LATER!!!)
  CX_L = 1010; % uatm,  Lower boundary for ingestion  (FOR DEVELPMENT ONLY, DELETE LATER!!!)
  CX_H = 6779; % uatm,  Upper boundary for ingestion  (FOR DEVELPMENT ONLY, DELETE LATER!!!)
  
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
  p_AmT = c_T * p_Am ;
  v_T   = c_T * v;  
  p_TT  = c_T * p_T;  
  k_JT  = c_T * k_J;
  kT_M  = c_T * k_M;
  p_XmT = p_AmT / kap_X;                  
  E_m   = p_AmT / v;                      % Jcm3, maximun reserve density
  
  % Test parameters (FOR DEVELPMENT ONLY, DELETE LATER!!!))
  E_m   = p_AmT / v;                  % Jcm3, maximun reserve density
  e     = (E / V) / E_m;              % -,  scaled energy density
  L_Vm =  v_T / (g * kT_M);           % cm, maximum volumetric lenght 
  
  CM_H = 1300;  % uatm, Higher pCO2 boundary for maintenace cost           (DEVELPMENT ONLY, DELETE LATER!!!)
  CM_L = 700;   % uatm, Lower pCO2 boundary for maintenacne cost           (DEVELPMENT ONLY, DELETE LATER!!!)
  q = 0.009;    % Jcm2 uatm-1 d-1, maintenace coefficient related to pCO2  (DEVELPMENT ONLY, DELETE LATER!!!)
  V_l = simu.V_l;

  %% Fluxes for the juvenile
  if E_H < E_Hb
      pX = 0;  % embryo stage
  else
      pX = f * p_XmT * V^(2/3);        % feeding 
  end
  pA = kap_X * pX * c_X;                     % assimilation (eq. 2.3 p.35)
  pM = (p_M * V + q * min(max(CO2 - CM_L,0), CM_H - CM_L * V^(2/3))) * c_T;  % Somatic maintenance 
  pJ = min(V, V_l) * ((1 - kap) / kap) * p_M * c_T;                          % Maturiry maintenance (eq.2.19 p.47)
  pC = (E/V) * ((E_G * v * V^(2/3) + p_M * V) / (kap * E/V + E_G )) * c_T;   % Reserve mobilization
  pG = kap * pC - pM;                  % growth
  pR = (1 - kap) * pC - pJ;             % reproduction (eq.2.55 p.67)
  pD = pM + pJ + (1 + kap_R) * pR;      % dissipation
  
  % coefficients for pA, pD and pG
  vA = 0;
  vD = 2.37E-04;
  vG = 0.2;
  aG = 1;
  aD = 0;
 
  %% Differential equations
  dE = pA - pC;                                       % dE/dt energy reserves
  dV = pG / E_G;                            % dV/dt structural volume 
  dR = pR;                                  % resproduction buffer
                             
  if e < (V^(1/3) / L_Vm)                             % starvation condition
    dV = 0;  % dV/dt structural volume
    dR = 0;
  end
  
  dS = (1 / p_Am) * (vG * pG + vD * pD);               % dS/dt shell volume  
  dO = (aG * vG * pG + aD * vD *pD) / (vG*pG + vD*pD); % dO/dt shell opacity (calcified structure)
  
  if E_H < E_Hp
      dH = (1 - kap) * pC - pJ; % dEH/dt
      dR = 0; % dER/dt
  else
      dH = 0;
      dR = pR; % reproduction
  end
  
  d = [dE; dV; dH; dR; dS; dO]; 