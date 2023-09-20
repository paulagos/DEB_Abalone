function obs = get_obs(simu)
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
  filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/results_Haliotis_tuberculata2.mat';
  res = load(filename);
  par  = res.par;          % get primary parameters
  cPar = parscomp_st(par); % get compound parameters
  vars_pull(par);
  vars_pull(cPar);

%% State variables
t   = simu.tEVHRSO(:,1); % days,     Time
E   = simu.tEVHRSO(:,2); % Joules,   Reserve energy
V   = simu.tEVHRSO(:,3); % cm^3,     Structural volume
E_R = simu.tEVHRSO(:,5); % Joules,   Reproduction buffer
V_S = simu.tEVHRSO(:,6); % cm^3,     Shell structural volume

% Non-observable quantities
E_m = cPar.p_Am ./ par.v;                % [Em], J m^-3, maximun reserve density
e   = E ./ (V * E_m);                    % e,          -, scaled reserve density

% Observable quantities (lenghts)
L_w       = V.^(1/3) / del_M;           % cm, physical length

% Observable quantities (total wet-weights)
Ww_V      =  d_V * V + (E + E_R) * w_E / mu_E;  % g, total wet-weight

% Observable quantities (Energy)
E_V       = mu_V * d_V / w_V * V;        % Joules, Energy content of structure
E_w       = (E_V + E + E_R) ./ Ww_V;     % J/g, Total energy content 

% Observable quantities (eggs number)
F         = kap_R * E_R / simu.EVHRSO_init(1);  % Fecundity = egg number = kap_R * E_R / E_0

% Observable quantities (Gonad wet weight)
G_Ww      =  E_R * (w_E / (d_E * mu_E));          % g, gonad wet-weight

obs = [t, e, L_w, Ww_V, E_w, F, G_Ww];