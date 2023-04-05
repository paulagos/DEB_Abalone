function tEVHRSO = indiv(simu)
  %---------------------------------------------------------------
  % Compute model predictions (numerical integration)
  %  
  % simu: 1-structure with individual features (parameters, env...),;
  %             see init.m
  %
  % tEVHRS : (nt,4) matrix with time and state variables
  %         t, d,      time
  %         E , J,     reserve energy
  %         V , cm^3,  structural volume
  %         E_H , J ,  cumulated energy inversted into maturity
  %         E_R , J,   reproduction buffer
  %         E_S, cm^3, shell structural volume
  %         O  ,-,     shell opacity
  %         
  % called by: main.m
  % calls    :flux.m
  %
  % 2013/03/15 - Laure Pecquerie
  % 2023/02/22 - Paulo F. Lagos
  %--------------------------------------------------------------
  
  %% unpack parameters values
  filename = '/Users/lagosjime/Documents/Work/UBO - Brest/Haliotis_tuberculata/results_Haliotis_tuberculata.mat';
  res = load(filename);
  par  = res.par;          % get primary parameters
  cPar = parscomp_st(par); % get compound parameters
  vars_pull(par); 
  vars_pull(cPar);
 
  %%
  tc = simu.t0; % set current time to initial time
  EVHRSO_tc = simu.EVHRSO_init; % set initial values of state variables for the egg 
  tEVHRSO = zeros(0,7); % initialize outputs
      
  i = 0; % initialize year index
  while tc < simu.tm
      i = i + 1; % spawning every year
      tnext = simu.ts + 365 * i; % integrate between two spawning dates (one spawning each year)
      if tnext > simu.tm
          tnext = simu.tm;
      end
      
      t = unique([tc, ceil(tc): ceil(tnext)])'; 
      [t, EVHRSO] = ode23s(@(t, EVHRSO) flux(t, EVHRSO, simu, par, cPar), t, EVHRSO_tc); % integrate
      tEVHRSO = [tEVHRSO; [t, EVHRSO]]; % append to output
      E_Hc = EVHRSO(end, 3);
      
      if E_Hc >= E_Hp %; E_Hp
          EVHRSO(end,4) = 0;
      end
      
      EVHRSO_tc = EVHRSO(end,:)'; 
      tc = t(end);
      
  end   