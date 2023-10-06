function simu = init(RCP, varEnv, envFoodQ, gen_number)
  %---------------------------------------------------------------
  % Define simulation and individual features
  %  
  % simu: structure with individual features (parameters, env...)
  %        (simu(i).par = parameters of individual i)
  %
  % called by : main.m
  % calls : set_par.m
  %
  % 2013/03/15 - Laure Pecquerie
  % 2023/02/01 - Paulo F. Lagos
  %--------------------------------------------------------------
  
  %% Load and unpack parameters
  filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';
  res = load(filename);
  par  = res.par;          % get primary parameters
  cPar = parscomp_st(par); % get compound parameters
  vars_pull(par); 
  vars_pull(cPar);
  
  %% set inital simulation times
  t0 = 0;             % initial time % January 1st
  tm = 20 * 365;      % final time 
  ts = t0;            % spawning time 
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    Choose if environment is constant or variable and RCP scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  envT    = varEnv;   % 1 = T constant,    2 = T variable
  envX    = varEnv;   % 1 = constant food, 2 = variable food
  envpCO2 = varEnv;   % 1 = contant pCO2,  2 = variable pCO2
  envRCP  = RCP;
  
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Color for plots
 if     envT == 1 && envX == 1 && envpCO2 == 1
	col = 'k'; 
 elseif envT == 2 && envX == 2 && envpCO2 == 2
	col = 'r';
 elseif envT == 1 && envX == 2 && envpCO2 == 1
    col = 'g';
 elseif envT == 2 && envX == 2 && envpCO2 == 1
    col = 'c';
 elseif envT == 2 && envX == 1 && envpCO2 == 1
    col = 'm';
 else
    col = 'b';
 end

  %% initial forcing variable
  T = 15 + 273.15;
  f = 0.9;
  pCO2 = 420;
  
  % Initial states
  pars_UE0 = [V_Hb; g; k_J; k_M; v];
  [U_E0, ~, ~] = initial_scaled_reserve(f, pars_UE0);
  
  E_0  = U_E0 * p_Am;            % J, initial reserve energy E_0 = U_E0 * p_Am 
  V_0  = 1e-8;                   % cm^3, initial structural volume
  E_H0 = 0;                      % J, initial cum. energy invested into maturation
  E_R0 = 0;                      % J, initial reproduction buffer
  K    = p_Am / (kap_X * F_m);   % same unit as food, half-saturation coefficient
 
  %simu.par = set_par;
  simu.t0      = t0;
  simu.tm      = tm;
  simu.ts      = ts;
  simu.envT    = envT;
  simu.RCP     = envRCP;
  simu.foodQ   = envFoodQ;
  simu.envX    = envX;
  simu.envpCO2  = envpCO2;
  simu.Tinit    = T;
  simu.Xinit    = f.* K ./ (1 - f);
  simu.pCO2init  = pCO2;
  simu.EVHR_init = [E_0 ; V_0 ; E_H0; E_R0];
  simu.tEVHR     = zeros(tm+1,length(simu.EVHR_init)+1, gen_number);
  simu.col = col;
