function simu = init
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
 
  t0 = 1;             % initial time % January 1st
  tm = 15 * 365;      % final time  
  ts = t0;            % spawning time
   
  simu.par = set_par; % load parameters for simulation
  
  % Choose if environment is constant or variable
  
  envT   = 2;   % 1 = T constant,    2 = T variable
  envX   = 2;   % 1 = constant food, 2 = variable food
  envCO2 = 2;   % 1 = contant pCO2,  2 = variable pCO2
  
  % Color for plots
 if     envT == 1 && envX == 1 && envCO2 == 1
	col = 'k'; 
 elseif envT == 2 && envX == 2 && envCO2 == 2
	col = 'r' ;
 elseif envT == 1 && envX == 2 && envCO2 == 1
    col = 'g';
 elseif envT == 2 && envX == 2 && envCO2 == 1
    col = 'c';
 elseif envT == 2 && envX == 1 && envCO2 == 1
    col = 'm';
 else
    col = 'b';
 end

  % initial forcing variable
  T = 20 + 273.15;
  f = 0.99;
  CO2 = 380;
  
  % Initial states
  pars_UE0 = [simu.par(23),simu.par(24), simu.par(25),simu.par(26), simu.par(6)]; %[V_Hb; g; k_J; k_M; v];
  [U_E0, Lb, info] = initial_scaled_reserve(f, pars_UE0);
  
  L_b = Lb;
  X = ['Lenght at birth (cm) = ', num2str(L_b), num2str(info)]; % Leght at birth
  disp(X); 

  
  E_0  = U_E0 * simu.par(5);     % J, initial reserve energy E_0 = U_E0 * p_Am 
  V_0  = 1e-8;                   % cm^3, initial structural volume
  E_H0 = 0;                      % J, initial cum. energy invested into maturation
  E_R0 = 0;                      % J, initial reproduction buffer
  E_S0 = 0;                      % cm^3, initial shell structural volume
  O_P0 = 0;                      % initial shell opacity
 
  %simu.par = set_par;
  simu.t0 = t0;
  simu.tm = tm;
  simu.ts = ts;
  simu.envT    = envT;
  simu.envX    = envX;
  simu.envCO2  = envCO2;
  simu.Tinit   = T;
  simu.Xinit   = f.* simu.par(22) ./ (1 - f);
  simu.CO2init = CO2;
  simu.EVHRSO_init = [E_0 ; V_0 ; E_H0; E_R0; E_S0; O_P0];
  simu.col = col;
