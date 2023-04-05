function pars = set_par
  %--------------------------------------------------------------- 
  % Parameters (values, unit and description)
  % pars: n-vector of parameters
  %    
  % called by : init.m
  %
  % 2013/03/15 - Laure Pecquerie
  % 2023/02/02 - Paulo F. Lagos
  % notation corresponds to the add_my_pet routines
  % Reference: <http://www.bio.vu.nl/thb/deb/deblab/add_my_pet/add_my_pet.pdf *add_my_pet manual*> 

  %---------------------------------------------------------------
  filename = '/Users/lagosjime/Documents/Work/UBO - Brest/Haliotis_tuberculata/results_Haliotis_tuberculata.mat';
  res = load(filename);
  par  = res.par;          % get primary parameters
  cPar = parscomp_st(par); % get compound parameters
 
  % temperature correction
  T1 = par.T_ref;        % 1, K, Reference temperature ; 
  TA = par.T_A;          % 2, K, Arrhenius temperature ;
   
  % feeding & assimilation
  F_m   = par.F_m;       % 3, l/d.cm^2, {F_m} max spec searching rate
  kap_X = par.kap_X;     % 4, -, digestion efficiency of food to reserve
  p_Am  = cPar.p_Am;     % 5, J/cm^2/d, maximum surface-specific assimilation rate

  % mobilisation, maintenance, growth & reproduction
  v     = par.v;         % 6, cm/d, energy conductance
  kap   = par.kap;       % 7, -, allocation fraction to soma = growth + somatic maintenance
  kap_R = par.kap_R;     % 8, -, reproduction efficiency
  p_M   = par.p_M;       % 9, J/d.cm^3, [p_M], vol-specific somatic maintenance
  p_T   = par.p_T;       % 10, J/d.cm^2, {p_T}, surface-specific som maintenance
  k_J   = par.k_J;       % 11, 1/d, maturity maint rate coefficient
  E_G   = par.E_G;       % 12, J/cm^3, [E_G], spec cost for structure

  % life stages: E_H is the cumulated energy from reserve invested in maturation
  E_Hb = par.E_Hb;      % 13, J, E_H^b, maturity at birth
  E_Hp = par.E_Hp;      % 14, J, E_H^p, maturity at puberty

  % param to compute observable quantities
  del_M = par.del_M;    % 15, -, adults and juveniles shape coefficient (convert vol.length to physical length)
  d_V   = par.d_V; 	    % 16, g/cm^3, specific density of structure (dry weight)
  mu_V  = par.mu_V;     % 17, J/mol
  mu_E  = par.mu_E;     % 18, J/mol
  w_V   = cPar.w_V;     % 19, g/mol
  w_E   = cPar.w_E;     % 20, g/mol
  w     = cPar.w;       % 21  wet / dry weight coefficient
  
  % Additional parameters
  K     = cPar.p_Am / (par.kap_X * par.F_m); % 22, same unit as food, half-saturation coefficient
  V_Hb  = cPar.V_Hb;     % 23         
  g     = cPar.g;        % 24
  k_J   = par.k_J;       % 25
  k_M   = cPar.k_M;      % 26
  del_L = par.del_L;     % 27, -, larvae shape coefficient 
  E_Hj  = par.E_Hj;      % 28, J, maturity at metamorphosis 
   
  %% pack parameters
  pars = [T1;   TA;   F_m; kap_X; p_Am; v;    kap;    kap_R; ...
          p_M;  p_T;  k_J; E_G;   E_Hb; E_Hp; del_M;  d_V; ...
          mu_V; mu_E; w_V; w_E;   w; K; V_Hb; g; k_J; k_M; del_L;...
          E_Hj];
