function GSI =  tGSI(t)
%---------------------------------------------------------------
  % Determine gonado somatic index  at time t, and temperature T
  % Here we compute it for a one year cycle temperature 
  % in the Bay of Biscay ())
  %  
  % t : scalar or n-vector of time points
  %
  % X : scalar or n-vector with food density (proxy = primary production)
  %
  % called by : main.m
  % calls : fnfourier.m (from DEBtool)
  % 
  % 2023/03/13 - Paulo F. Lagos
  %--------------------------------------------------------------
tTGSI = [...% time (d) ~ Temperature (deg C) ~ gonado somatic index (GSI)
30.4	8.89	41.26174795
121.7	10.35	54.27255236
182.6	13.85	74.15858803
185.6	13.96	85.66824539
213.1	15.78	77.32257826
243.5	18.56	72.53311626
273.9	16.87	40.1446755
334.8	13.91	39.54044254
365.2	11.92	32.87463924
];
tGSI = tTGSI(:,[1 3]);   % time ~ GSI data
tT   = tTGSI(:,[1 2]);   % time ~ temperatures for one year period data
  
  % time - GSI: max gonad weight and proportion to soma
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
 [~, ~, tau_b l_j l_p l_b l_i , ~, ~, info] = get_tj(pars_tj, f);
  L_b = L_m * l_b;
  L_j = L_m * l_j;
  E_b = f * L_b^3 * p_Am / v ;
  EaT_b = tau_b / k_M / tempcorr(C2K(13.7), T_ref, T_A); % average temperature for North of Brittany 
  ELH_0 = [E_b; L_b; E_Hb; 0];                           % initial conditions for state variables E,L,H and R
  
  tGSI(:,2) = C2K(tGSI(:,2)); % convert to K before the ode-solver for efficiency
  [~, ELHR] = ode45(@dget_ELHR_j, [EaT_b; tGSI(:,1)], ELH_0, [], L_b, L_j, L_m, p_Am, v, g, k_J, kap, kap_R, E_Hb, E_Hj, E_Hp, T_ref, T_A, tGSI, f); %[] in third column means no ODE options
  ELHR(1,:)=[];
  E = ELHR(:,1);   % reserve E,  J
  L = ELHR(:,2);   % structural length, cm
  R = ELHR(:,4);   % reproduction buffer E_R, J
  % H = ELHR(:,3); % maturity E_H, J
  W = L.^3 + (E + R) * w_E / mu_E / d_E;    % somatic weight (no shell)
  GSI = (R * w_E / mu_E / d_E) ./ W .* 100; % GSI