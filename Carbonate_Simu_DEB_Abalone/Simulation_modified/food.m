function X = food(t)
%---------------------------------------------------------------
  % Determine food density at time t,
  % Here we compute it with a fourier transform of a primary production cycle
  % in the Bay of Biscay (Pecquerie et al. 2009)
  %  
  % t : scalar or n-vector of time points
  %
  % X : scalar or n-vector with food density (proxy = primary production)
  %
  % called by : main.m
  % calls : fnfourier.m (from DEBtool)
  % 
  % 2013/03/15 - Laure Pecquerie
  %--------------------------------------------------------------
   
 p5 = [365.0000  146.2594;
       53.5370  -12.0523;
      -46.2953  -22.1553;
        6.5684   16.8142;
        13.5536    4.0526;
        -0.6571   -7.5133];

X = fnfourier(t, p5);