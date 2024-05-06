function X = food(t,simu)
%---------------------------------------------------------------
  % Determine food density at time t,
  % computed it with a fourier transform of anual cycle of
  % total algae biomass (kg m2 h-1) in the area of Aber Wrac'h
  %  
  % t : scalar or n-vector of time points
  %
  % X : scalar or n-vector with food density (proxy = primary production)
  %
  % called by : main.m
  % calls : fnfourier.m (from DEBtool)
  % 
  % 2024/02/05 - Paulo F. Lagos
  %--------------------------------------------------------------
%   
%   p5 = [0   0.120336841
%         25	0.120336943
%         51	0.080808081
%         61	0.121739945
%         62	0.153388937
%         64	0.045127808
%         79	0.128389752
%         80	0.133067199
%         81	0.171356421
%         92	0.143891362
%         102	0.213522427
%         107	0.114707087
%         111	0.082976792
%         112	0.066100680
%         118	0.516795866
%         121	0.051078320
%         125	0.083091169
%         126	0.096259626
%         127	0.283224401
%         128	0.016692169
%         129	0.007280964
%         134	0.005624692
%         136	0.229242782
%         137	0.124843945
%         139	0.305610129
%         140	0.179349161
%         141	0.009272997
%         145	0.239533997
%         148	0.018570923
%         149	0.012698413
%         154	0.011643035
%         155	0.151752264
%         158	0.147601476
%         159	0.014402881
%         164	0.006878630
%         167	0.134374738
%         168	0.226372383
%         174	0.063268893
%         177	0.065956367
%         178	0.096982721
%         193	0.026360384
%         197	0.035219018
%         204	0.061957869
%         207	0.073003358
%         217	0.098645861
%         228	0.143860540
%         234	0.181950509
%         235	0.294117647
%         236	0.198603441
%         237	0.257400257
%         241	0.127029495
%         242	0.107642626
%         243	0.02085673
%         251	0.020030817
%         252	0.020701076
%         253	0.028871391
%         254	0.013310576
%         259	0.011726079
%         260	0.013291981
%         261	0.039114398
%         263	0.022681452
%         269	0.007715454
%         272	0.003639206
%         279	0.017107205
%         289	0.003710667
%         290	0.013002492
%         295	0.002443176
%         301	0.025974026
%         311	0.007279465
%         327	0.021600000
%         329	0.130718954
%         339	0.063356300
%         340	0.014102715
%         349	0.013049625
%         356	0.001113117
%         357	0.003595075
%         358	0.011579100
%         365	0.005972526];


%My coefficients from algal data
p5 = [365.000   0.0861;
      0.0205   -0.0471;
     -0.0082    0.0175;
     -0.0177    0.0203;
     -0.0274    0.0088;
     -0.0010   -0.0293];
     
% t1 = p5(:,1);    
% t = linspace(0, 365,100)'; % time points for evaluation
% fr = 0.1;% function value to solve x from f(x) = fr

%% spline1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [fs1, dfs1] = spline1(t1, p5); % function values and derivative
% ifs1 = ispline1(t1, p5); % integration
% ts1_fr = rspline1(p5,0,0,fr); % roots of f(x) = fr;
% [ts1_min, ts1_max] = espline1(p5); % extremes


%% spline, here natural spline. 

% [fs, dfs, ddfs] = spline(t1, p5); % function values and derivatives
% ifs = ispline(t1, p5); % integration
% ts_fr = 0; ts_min = 0; ts_max = 0;
%     
% % get Fourier coefficients with 5 terms, s0 1 + 2 * 4 = 9 pars
% p = get_fourier(365, 5, p5); % 365 = period(days), 
% ff = fnfourier(t1, p);  % function values
% dff = dfnfourier(t1,p); % derivative
% ddff = ddfnfourier(t1,p); % second derivative
% iff = ifnfourier(t1,p);  % integral
% tf_fr = 0; tf_min = 0; tf_max = 0;

%
X = fnfourier(t + simu.tb , p5);
% 
X = max(X, 0.05); % Make sure that we dont have negative food values coming  from the fourier transform
% I COULD USE ANOTER MINUM VALUE THAT IS NOT 0 IF i BELOEVE THE ABALONE
% NEVER REACH STARVATION

% % Food in the environment plots 
% subplot(2,3,1)
% plot(p5(:,1), p5(:,2), '*k', t1, fs1, 'b', ...
%     [0;365], [fr;fr], 'm', ts1_fr, fr + 0 * ts1_fr, 'om', ...
%     ts1_min(:,1), ts1_min(:,2), '+b', ts1_max(:,1), ts1_max(:,2),'+r')
% xlabel('independent var (time)')
% ylabel('dependent var')
% title('first order spline')
% 
% subplot(2,3,2)
% plot(p5(:,1), p5(:,2), '*k', t1, fs, 'r', ...
%     [0;365], [fr;fr], 'm', ts_fr, fr + 0 * ts_fr, 'om', ...
%     ts_min(:,1), ts_min(:,1), '+b', ts_max(:,1), ts_max(:,1),'+r')
% xlabel('independent var')
% ylabel('dependent var')
% title('cubic spline')
% 
% subplot(2,3,3)
% plot(p5(:,1), p5(:,2), '*k', t1, ff, 'g', ...
%     [0;365], [fr;fr], 'm', tf_fr, fr + 0 * tf_fr, 'om', ...
%     tf_min(:,1), tf_min(:,1), '+b', tf_max(:,1), tf_max(:,1),'+r')
% xlabel('independent var')
% ylabel('dependent var')
% title('Fourier series')
% 
% 
% subplot(2,3,4)
% plot(t, ifs1, 'b', t1, ifs, 'r', t1, iff,'g')
% xlabel('independent var')
% ylabel('integrated')
% 
% subplot(2,3,5)
% plot(t, dfs1, 'b', t1, dfs, 'r', t1, dff, 'g')
% xlabel('independent var')
% ylabel('derivative')
% 
% subplot(2,3,6)
% plot([0; 365], [0;0], 'b', t1, ddfs, 'r', t1, ddff, 'g')
% xlabel('independent var')
% ylabel('second deribative')

end
 