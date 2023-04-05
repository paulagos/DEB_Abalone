function obs = get_obs(simu)
  %---------------------------------------------------------------
  % Compute physical length, weight, fecundity, and energy content from state variables
  %  
  % obs : (n, 5) matrix with 
  %         t: d, time,
  %			Lw: cm, physical length
  %         Ww: g, total wet weight,
  %			F: #, fecundity
  %         Ew: J/g, energy content per unit wet weight
  %     
  %
  % called by : main.m
  % 
  % 2013/03/15 - Laure Pecquerie
  % 2023/03/21 - PaulI F. Lagos
  %-------------------------------------------------------------
%% unpack parameters, state variables
 filename = '/Users/lagosjime/Documents/Work/UBO - Brest/Haliotis_tuberculata/results_Haliotis_tuberculata.mat';
  res = load(filename);
  par  = res.par;          % get primary parameters
  cPar = parscomp_st(par); % get compound parameters
  vars_pull(par); 
  vars_pull(cPar);

%% State variables
t   = simu.tEVHRSO(:,1); % Time
E   = simu.tEVHRSO(:,2); % Reserve
V   = simu.tEVHRSO(:,3); % Structural volume
E_R = simu.tEVHRSO(:,5); % Reproductive buffer
V_s = simu.tEVHRSO(:,6); % Shell volume

% Observable quantities 
L_w       = V.^(1/3) / del_M;           % cm, physical length
W_V       = d_V * V;                    % dry weight of structure
W_E_and_R = w_E / mu_E * (E + E_R) ;    % dry weight of E and E_R
W         = W_V + W_E_and_R;            % total dry weight
W_w       = w * W;                      % assume the same water content in structure and reserves
E_V       = mu_V * d_V / w_V * V ;  
E_w       = (E_V + E + E_R) ./ W_w ;
F         = kap_R * E_R / simu.EVHRSO_init(1);  % Fecundity = egg number = kap_R * E_R / E_0

% Observable quantities  for the shell
S_l       = V_s.^(1/3) / del_M;         % cm, Shell physical lenght

obs = [t , L_w , W_w, E_w, F, S_l];