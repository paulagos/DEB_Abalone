function obs = get_obs(simu, zoom, gen)
  %---------------------------------------------------------------
  % Compute physical length, weight, fecundity, and energy content from state variables
  %  
  % obs : (n, 8) matrix with 
  %                 (1) = t,   d,  time , 
  %                 (2) = L_w, cm, physical length, 
  %                 (3) = S_l, cm, Shell physical lenght 
  %                 (4) = W_S, g,  dry weight of the shell, 
  %                 (5) = W,   g,  total dry weight (tissue +  reserves + gonad + shell)
  %                 (6) = W_w, g,  total wet weight, assume the same water content in structure and reserves 
  %                 (7) = E_w, J/g,total energy content 
  %                 (8) = F,   #,  Fecundity or total number of eggs 
  %  
  % called by : main.m
  % 
  % 2013/03/15 - Laure Pecquerie
  % 2023/03/21 - PaulI F. Lagos
  %-------------------------------------------------------------
  
%% State variables
t   = simu.tEVHR(:,1, zoom, gen); % days,     Time
E   = simu.tEVHR(:,2, zoom, gen); % Joules,   Reserve energy
V   = simu.tEVHR(:,3, zoom, gen); % cm^3,     Structural volume
E_R = simu.tEVHR(:,5, zoom, gen); % Joules,   Reproduction buffer

% Observable quantities
del_M     =  286.2^(1/3) / 11.25; % ultimate total lenght and weight using data from my_data file Page 11 DEB book
L_w       =  V.^(1/3) / del_M; %simu.par.del_M;                  % cm, physical length
% L_w       =  V.^(1/3) / simu.par.del_M;                  % cm, physical length
Ww_V      =  simu.par.d_V * V + (E + E_R) * simu.cPar.w_E / simu.par.mu_E;  % g, total wet-weight
F         =  simu.par.kap_R * E_R / simu.EVHR_init(1);   % Fecundity = egg number = kap_R * E_R / E_0
G_Ww      =  E_R * simu.cPar.w_E / simu.par.mu_E;  % g, gonad wet-weight
GSI       =  G_Ww ./ (Ww_V - G_Ww);              % -, Gonado somatic index (GSI)

% Non-observable quantities
E_m       = simu.cPar.p_Am ./ simu.par.v;          % [Em], J m^-3, maximun reserve density
e         = E ./ (V * E_m);              % e,          -, scaled reserve density
E_V       = simu.par.mu_V * simu.par.d_V / simu.cPar.w_V * V;        % Joules, Energy content of structure
E_w       = (E_V + E + E_R) ./ Ww_V;     % J/g, Total energy content 

obs = [t, e, L_w, Ww_V, E_w, F, G_Ww, GSI];