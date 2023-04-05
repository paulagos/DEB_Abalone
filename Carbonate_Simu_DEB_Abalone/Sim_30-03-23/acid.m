function CO2 = acid(t)
%---------------------------------------------------------------
  % Determine pCO2 at time t, 
  % Here we compute it as a sinus function that is ajusted
  % to a pH cycle (water column) in the Bay of Biscay 
  %  
  % t : scalar or n-vector of time points
  %
  % pCO2 : scalar or n-vector with partial pressure of pCO2
  %
  % called by : main.m
  % 
  % 2013/03/15 - Paulo F. Lagos
  %--------------------------------------------------------------
  
  CO2 = 380 + 80 * sin(2 * pi * (t + 107) / 365);
