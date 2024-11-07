function simu = init(filename, SSP, location, varEnv, envFoodQ, gen_number, zoom, born_date) %, pco2, temp)
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
  res = load(filename); 
  par  = res.par;          % get primary parameters
  cPar = parscomp_st(par); % get compound parameters

  % Set parameter variables
  [V_Hb, g, k_J, k_M, v]    = deal(cPar.V_Hb, cPar.g, par.k_J, cPar.k_M, par.v);
  [p_Am] = deal(cPar.p_Am);
  
  % set simulation times
  t0 = 0;     % initial time % 0; January 1st
  tm = 20 * 365;      % max simulation length in number of days
  ts = t0;            % spawning time xs
  
  % Set environment flags
  envT    = varEnv;  
  envX    = varEnv;   
  envpCO2 = varEnv;   
  envSSP  = SSP;
  
  % Define plot color based on environment
  col = getColor(envT, envX, envpCO2);
  
  % constant environemantal conditions
  T    = 14 + 273.15; %14 + 273.15 or  %temp + 273.15; %
  f    = 0.9;
  pCO2 = 420; % pco2 ;
  
  % SSP scenarios applied to in-situ environmental data
  [ssp_temp, ssp_pco2] = forcing_deb(location, SSP, born_date);
   
  % Calculate initial states
  pars_UE0     = [V_Hb; g; k_J; k_M; v];
  [U_E0, ~, ~] = initial_scaled_reserve(f, pars_UE0);
  E_0          = U_E0 * p_Am;            % J, initial reserve energy E_0 = U_E0 * p_Am
  V_0          = 1e-8;                   % cm^3, initial structural volume
  E_H0         = 0;                      % J, initial cum. energy invested into maturation
  E_R0         = 0;                      % J, initial reproduction buffer
  K            = 0.01;                   % p_Am / (kap_X * F_m)  same unit as food, half-saturation coefficient
        
  % Set simulation parameters in the 'simu' structure
  simu.t0        = t0;
  simu.tm        = tm;
  simu.ts        = ts;
  simu.tb        = born_date;
  simu.envT      = envT;
  simu.SSP       = envSSP;
  simu.foodQ     = envFoodQ;
  simu.envX      = envX;
  simu.envpCO2   = envpCO2;
  simu.zoom      = zoom;
  simu.par       = par;
  simu.cPar      = cPar;
  simu.Tinit     = T;
  simu.Xinit     = f.* K ./ (1 - f);
  simu.pCO2init  = pCO2;
  simu.K         = K;
  simu.tempData  = ssp_temp;
  simu.pco2Data  = ssp_pco2;
  simu.EVHR_init = [E_0 ; V_0 ; E_H0; E_R0];
  simu.tEVHR     = zeros(tm+1,length(simu.EVHR_init)+1, gen_number);
  simu.col       = col;
  
end
  
  function col = getColor(envT, envX, envpCO2)
    % Determine plot color based on environment conditions
    if envT == 1 && envX == 1 && envpCO2 == 1
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
    
  end
