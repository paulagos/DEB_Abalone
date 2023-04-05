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
  

  par = simu.par;
  tc = simu.t0; % set current time to initial time
  EVHRSO_tc = simu.EVHRSO_init; % set initial values of state variables for the egg 
  
  tEVHRSO = zeros(0,7); % initialize outputs
  
  odeoptions = odeset('Event', @eventfunc_deb);
  
  %t = [tc:(1/24):10]'; 
  t = [tc:10]';
      [t, EVHRSO, TE, YE, ie] = ode23s(@flux_egg, t, EVHRSO_tc, odeoptions, simu);
      tEVHRSO = [tEVHRSO; [t, EVHRSO]]; % append to output
      if isempty(ie)
          print('the egg is not hatching, try again')
          return
      end
      tc = t(end);
      EVHRSO_tc = EVHRSO(end,:)'; % set initial values of state variables for the larvae
      V_egg = EVHRSO(end,2);
      simu.V_egg = V_egg;         % gives final value of state variable V for the egg
      
      
      
 t = [tc, ceil(tc):(ceil(tc) + 50)]'; 
      [t, EVHRSO, TE, YE, ie] = ode23s(@flux_larvae, t, EVHRSO_tc, odeoptions, simu);
      tEVHRSO = [tEVHRSO; [t, EVHRSO]]; % append to output
      tc = t(end);
      if isempty(ie)
          disp('the larvae is not transitioning to juvenile, try again')
          return
      end
      EVHRSO_tc = EVHRSO(end,:)'; % set initial values of state variables for the larvae
      V_l = EVHRSO(end,2);
      simu.V_l = V_l;              % gives final value of state variable V for the larvae
      
 t = [tc, ceil(tc): (ceil(tc) + 2150)]'; 
      [t, EVHRSO, TE, YE, ie] = ode23s(@flux_juv, t, EVHRSO_tc, odeoptions, simu);
      tEVHRSO = [tEVHRSO; [t, EVHRSO]]; % append to output 
      tc = t(end);
      if isempty(ie)
          disp('the juvenile is not transitioning to adult, try again')
          return
      end
      EVHRSO_tc = EVHRSO(end,:)';
      t_p = tc;
      V_p = EVHRSO(end,2);
      simu.V_p = V_p; 
  
  i = 0; % initialize year index
  while tc < simu.tm
      i = i + 1; % spawning every year
      tnext = t_p + simu.ts + 365 * i; % integrate between two spawning dates (one spawning each year)
      if tnext > simu.tm
          tnext = simu.tm;
      end
      
      t = unique([tc, ceil(tc): ceil(tnext)])'; 
      [t, EVHRSO] = ode23s(@(t, EVHRSO)flux_adult(t, EVHRSO, simu), t, EVHRSO_tc); % integrate
      tEVHRSO = [tEVHRSO; [t, EVHRSO]]; % append to output
      E_Hc = EVHRSO(end, 3);
      if E_Hc >= par(14) %; E_Hp
          EVHRSO(end,4) = 0;
      end
      
      EVHRSO_tc = EVHRSO(end,:)'; 
      tc = t(end);
      
  end   