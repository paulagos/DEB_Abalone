function d = flux_egg(t, EVHRSO, simu)
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
   
  %% unpack state vars
  E    = EVHRSO(1);  % J, reserve energy
  V    = EVHRSO(2);  % cm^3, structural volume
  E_H  = EVHRSO(3);  % J , cumulated energy inversted into maturity 
  E_R  = EVHRSO(4);  % J, reproduction buffer
  E_S  = EVHRSO(5);  % Shell structural volume
  O    = EVHRSO(6);  % shell opacity
  
  
  % read parameter values
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
  
  
  % simplest temperature correction function, 
  % see Kooijman 2010 for more detailed formulation (e.g. p. 21)
  c_T   = exp(TA/ T1 - TA/ T) ;
  p_AmT = c_T * p_Am ;
  v_T   = c_T * v;  
  p_MT  = c_T * p_M;                     %, somatic maintenance
  p_TT  = c_T * p_T;  
  k_JT  = c_T * k_J;
  kT_M  = c_T * k_M;      
  
  % Fluxes for the egg stage
  pX = 0;                              % feeding 
  pA = 0;                              % assimilation (eq. 2.3 p.35)
  pM = p_MT * V;                       % somatic mnaintenace 
  pJ = k_JT * E_H;                     % Maturiry maintenance (eq.2.19 p.47)
  pT = p_TT * V^(2/3);                 % Surface-related somatic maintenance
  pS = pM + pT;                        % Somatic maintenance (eq.2.17 p.42)
  pC = (E/V) * (E_G * v_T * V^(2/3) + pS ) / (kap * E/V + E_G ); %  Mobilisation (eq. 2.12 p.37)
  pG = kap * pC - pM;                  % growth
  pR = (1 - kap) * pC - pJ;             % maturation (eq.2.55 p.67)
  pD = 0;                               % dissipation
  
   
  % Differential equations
  dE = pA - pC;                                       % dE/dt energy reserves
  dV = (kap * pC - pS) / E_G;                         % dV/dt structural volume 
  dH = pR;                                            % maturity
  dR = 0;                                             % repro buffer
  dS = 0;                                             % dS/dt volumetric lenght 
  dO = 0;                                             % dO/dt shell opacity (calcified structure)
   
  d = [dE; dV; dH; dR; dS; dO]; 