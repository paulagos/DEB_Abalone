function obs = get_obs(filename, simu, gen)
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
  %% unpack parameters, state variables
  res = load(filename);
  par  = res.par;          % get primary parameters
  cPar = parscomp_st(par); % get compound parameters
  vars_pull(par);
  vars_pull(cPar);
  
%% State variables
t   = simu.tEVHR(:,1, gen); % days,     Time
E   = simu.tEVHR(:,2, gen); % Joules,   Reserve energy
V   = simu.tEVHR(:,3, gen); % cm^3,     Structural volume
E_R = simu.tEVHR(:,5, gen); % Joules,   Reproduction buffer

% Observable quantities
L_w       =  V.^(1/3) / del_M;                  % cm, physical length
Ww_V      =  d_V * V + (E + E_R) * w_E / mu_E;  % g, total wet-weight
F         =  kap_R * E_R / simu.EVHR_init(1);   % Fecundity = egg number = kap_R * E_R / E_0
G_Ww      =  E_R * (w_E / (d_E * mu_E));        % g, gonad wet-weight
GSI       =  G_Ww ./ (Ww_V - G_Ww);              % -, Gonado somatic index (GSI)

% Non-observable quantities
E_m       = cPar.p_Am ./ par.v;          % [Em], J m^-3, maximun reserve density
e         = E ./ (V * E_m);              % e,          -, scaled reserve density
E_V       = mu_V * d_V / w_V * V;        % Joules, Energy content of structure
E_w       = (E_V + E + E_R) ./ Ww_V;     % J/g, Total energy content 

obs = [t, e, L_w, Ww_V, E_w, F, G_Ww, GSI];