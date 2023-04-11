function get_plots(simu)

filename = '/Users/lagosjime/Documents/Work/UBO - Brest/Haliotis_tuberculata/results_Haliotis_tuberculata.mat';
  res = load(filename);
  par  = res.par;          % get primary parameters
  cPar = parscomp_st(par); % get compound parameters

  % time , environmental forcing
  t = simu.tEVHRSO(:,1);
  
  if simu.envX == 1 
          X = repmat(simu.Xinit, 1, length(t));
  else
      X = food(t);
  end
  
  if simu.envT == 1
          T = repmat(simu.Tinit, 1, length(t)) - 273; % in DegC
  else
          
          T = temp(t) - 273.15;  
  end
  
  if simu.envCO2 == 1
      CO2 = repmat(simu.CO2init, 1, length(t));
      
  else 
      CO2 = acid(t);
  end
  
  %% state variables
  E   = simu.tEVHRSO(:,2); % reserves
  V   = simu.tEVHRSO(:,3); % structural volume
  E_H = simu.tEVHRSO(:,4); % maturity
  E_R = simu.tEVHRSO(:,5); % reproduction buffer
  E_S = simu.tEVHRSO(:,6); % shell volume
  O   = simu.tEVHRSO(:,7); % shell opacity
  
  %% quantities
  K     = cPar.p_Am / (par.kap_X * par.F_m); % 22, same unit as food, half-saturation coefficient
  f = X ./ (X + K); %scaled functional response
  Em = cPar.p_Am / par.v;
  sc_res_dens   = E ./ (V * Em); % -, scaled reserve density e
  
  % spawning date indices
  i_sp = find(and((E_R == 0),(E_H >= par.E_Hp)));
  i_sp = i_sp - 1 ; % look at the preceding line with E_R value before spawning
  
  %% observable quantities
  Lw  = simu.obs(:,2);
  Ww  = simu.obs(:,3);
  Ew  = simu.obs(:,4);
  F   = simu.obs(:,5);
  S_l = simu.obs(:,6); % shelll length\
  
 
  % graph 1
  nfig = 1; % number of figures
  for i = 1:nfig
    figure(i)
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperUnits','in');
    set(gcf,'PaperPosition',[0 0 14 10]); 
    set(gca,'Fontsize', 12)
  end
  
  figure(1)
  subplot(3, 4, 1)
  hold on
  plot( t, T, simu.col, 'Linewidth',1.5)
  title('Temperature')
  xlabel('Time (d)')
  ylabel('T (°C)')
  axis([0 t(end) 0 25])

  subplot(3, 4, 2)
  hold on
  plot( t, X , simu.col, 'Linewidth',1.5)
  title('Food')
  xlabel('Time (d)')	
  ylabel( 'Food density (units up to you)')
  
  subplot
  subplot(3, 4, 3)
  hold on
  plot( t, CO2, simu.col, 'Linewidth',1.5)
  title('pCO2 levels')
  xlabel('Time (d)')	
  ylabel( 'pCO2')
  
  
  subplot(3, 4, 4)
  hold on
  plot( t, sc_res_dens, [simu.col '-'], 'Linewidth',1.5)
  axis([0 t(end) 0 1.7])
  title('Scaled reserve density')
  xlabel('Time (d)')
  ylabel('e (-)')
  
  subplot(3, 4, 5)
  hold on
  plot(t, E,  simu.col, 'Linewidth',1.5)
  title('Reserve energy')
  xlabel('Time (d)')
  ylabel('E (J)')
  
  subplot(3, 4, 6)
  hold on
  plot(t, V, simu.col, 'Linewidth',1.5)
  title('Structural volume')
  xlabel('Time (d)')
  ylabel('V (cm^3)')
  
  subplot(3, 4, 7)
  hold on
  plot(t, E_H, simu.col, 'Linewidth',1.5)
  title('Energy invested into maturation')
  xlabel('Time (d)')
  ylabel('E_H (J)')
  
  subplot(3, 4, 8)
  hold on
  plot(t, E_R, simu.col , 'Linewidth',1.5)
  title('Reproduction buffer (J)')
  xlabel('Time (d)')
  ylabel('E_R (J)')
  
  subplot(3, 4, 9)
  hold on
  plot(t, Lw, simu.col , 'Linewidth',1.5)
  title('Physical length (cm)')
  xlabel('Time (d)')
  ylabel(' Lw (cm)')
 
  subplot(3, 4, 10)
  hold on
  plot(t, Ww, simu.col , 'Linewidth',1.5)
  title('Total Weight')
  xlabel('Time (d)')
  ylabel('Ww (g)')
  
  subplot(3, 4, 11)
  hold on
  plot(t, Ew, simu.col , 'Linewidth', 1.5)
  title('Energy content (J/g)')
  xlabel('Time (d)')
  ylabel('<Ew> (J.g^{-1})')

  subplot(3, 4, 12)
  hold on
  plot(Lw(i_sp), F(i_sp), simu.col, 'Linewidth',1.5)
  title('Fecundity-Length')
  xlabel('Lw (cm)')
  ylabel('F (#)')
 
  saveas(gcf, 'fig_DEB_standard.jpg', 'jpg') 
  % other ways to save the graph
  % saveas(gcf, 'fig_DEB_standard.eps', 'psc2')
  % print -djpeg90 fig_DEB_standard.jpg 


% graph 2
  nfig = 1; % number of figures
  for i = 1:nfig
    figure(i)
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperUnits','in');
    set(gcf,'PaperPosition',[0 0 14 10]); 
    set(gca,'Fontsize', 14)
  end
  
  figure(2)
  subplot(1, 3, 1)
  hold on
  plot(t, S_l, simu.col, 'Linewidth', 1.5)
  xlabel('Time (d)')
  ylabel('Shell lenght')
  v = axis;
  text((v(2)-v(1)) * 0.07 + v(1), (v(4)-v(3)) *.9 + v(3), '(a)','Fontsize',12)
  %
  subplot(1, 3, 2)
  hold on
  plot(t, O, simu.col, 'Linewidth', 1.5)
  xlabel('Time (d)')
  ylabel('Shell opacity')
  v = axis;
  text((v(2)-v(1)) * 0.07 + v(1), (v(4)-v(3)) *.9 + v(3), '(b)','Fontsize',12)
  %
  subplot(1, 3, 3)
  hold on
  plot(S_l, O, simu.col, 'Linewidth', 1.5)
  xlabel('Shell lenght (cm)')
  ylabel('Shell Opacity (-)')
  v = axis;
  text((v(2)-v(1)) * 0.07 + v(1), (v(4)-v(3)) *.9 + v(3), '(c)','Fontsize',12)
  %
  
  % temperature correction factor plots 
  T = C2K(T);
  T_ref = C2K(16);
  T_A1  = 9000;
  T_A2  = 7000;
  T_A3  = 5000;
  
  T_AH = 94000;
  T_AL = 59000;
  T_H  = 295;
  T_L  = 280;
  
  T_AH2 = 92000;
  T_AL2 = 50000;
  T_H2  = 293;
  T_L2  = 278;
 
  c_T1 = exp(T_A1 / T_ref - T_A1 ./ T);
  c_T2 = exp(T_A2 / T_ref - T_A2 ./ T);
  c_T3 = exp(T_A3 / T_ref - T_A3 ./ T);
  
    s_L_ratio = (1 + exp(T_AL/ T_ref - T_AL/ T_L)) ./ ...
	            (1 + exp(T_AL ./ T   - T_AL/ T_L));
    
    s_H_ratio = (1 + exp(T_AH/ T_H - T_AH/ T_ref)) ./ ...
	            (1 + exp(T_AH/ T_H - T_AH ./ T  ));
            
    s_L_ratio2 = (1 + exp(T_AL2/ T_ref - T_AL2/ T_L2)) ./ ...
	             (1 + exp(T_AL2 ./ T   - T_AL2/ T_L2));  
             
    s_H_ratio2 = (1 + exp(T_AH2/ T_H2 - T_AH2/ T_ref)) ./ ...
	             (1 + exp(T_AH2/ T_H2 - T_AH2 ./ T  ));
            
    T_C1 = c_T1 .* ((T <= T_ref) .* s_L_ratio + (T > T_ref) .* s_H_ratio);
    T_C2 = c_T2 .* ((T <= T_ref) .* s_L_ratio + (T > T_ref) .* s_H_ratio);
    T_C3 = c_T3 .* ((T <= T_ref) .* s_L_ratio + (T > T_ref) .* s_H_ratio);
    
    T_C4 = c_T1 .* ((T <= T_ref) .* s_L_ratio2 + (T > T_ref) .* s_H_ratio2);
    T_C5 = c_T2 .* ((T <= T_ref) .* s_L_ratio2 + (T > T_ref) .* s_H_ratio2);
    T_C6 = c_T3 .* ((T <= T_ref) .* s_L_ratio2 + (T > T_ref) .* s_H_ratio2);
  
  figure (3)
  subplot(1, 3, 1)
  hold on
  plot(K2C(T), c_T1 , 'Color',[0.9, 0.3, 0.0], 'LineWidth', 2);
  hold on
  plot(K2C(T), c_T2, 'Color',[0.6, 0.2, 0.9], 'LineWidth', 2);
  hold on
  plot(K2C(T), c_T3, 'Color',[0.4, 0.4, 0.3], 'LineWidth', 2);
  xlabel('Temperature ({\circ}C)')
  ylabel('Correction factor (c_T)')
  title('Temperature corr. factor at 16{\circ} ')
  legend('T_A = 9000', 'T_A = 7000', 'T_A = 5000');
  %
  subplot(1, 3, 2)
  hold on
  plot(K2C(T), T_C1 , 'Color',[0.9, 0.3, 0.0], 'LineWidth', 2);
  hold on
  plot(K2C(T), T_C2 , 'Color',[0.6, 0.2, 0.9], 'LineWidth', 2);
  hold on
  plot(K2C(T), T_C3 , 'Color',[0.4, 0.4, 0.3], 'LineWidth', 2);
  xlabel('Temperature ({\circ}C)')
  ylabel('Correction factor (c_T)')
  title('Temperature corr. factor at 16{\circ} ')
  subtitle ('T_{AH}=94000, T_{AL}=59000, T_H=295, T_L=280');
  legend('T_A = 9000', 'T_A = 7000', 'T_A = 5000');
  %
  subplot(1, 3, 3)
  hold on
  plot(K2C(T), T_C4 , 'Color',[0.9, 0.3, 0.0], 'LineWidth', 2);
  hold on
  plot(K2C(T), T_C5 , 'Color',[0.6, 0.2, 0.9], 'LineWidth', 2);
  hold on
  plot(K2C(T), T_C6 , 'Color',[0.4, 0.4, 0.3], 'LineWidth', 2);
  xlabel('Temperature ({\circ}C)')
  ylabel('Correction factor (c_T)')
  title('Temperature corr. factor at 16{\circ} ')
  subtitle ('T_{AH}=92000, T_{AL}=50000, T_H=293, T_L=278');
  legend('T_A = 9000', 'T_A = 7000', 'T_A = 5000');
  hold off
  %
  
   % Temperature corrected parameters
  p_AmT = c_T1 * cPar.p_Am;
  p_XmT = p_AmT / par.kap_X;
  
 % acidification factor correction factor plots
  p_M   = par.p_M;       % 9, J/d.cm^3, [p_M], vol-specific somatic maintenance
  CM_H2 = 15750;   % uatm, Higher pCO2 boundary for maintenace cost 
  CM_L2 = 250;     % uatm, Lower pCO2 boundary for maintenacne cost 
  q     = 5; % Jcm2 uatm-1 d-1, maintenace coefficient related to pCO2
   
  x = min(CO2):max(CO2);
  y = linspace(min(T), max(T),length(x));
  Z = zeros(length(x));
  
  c_T = exp(T_A1 / T_ref - T_A1 ./ y);
 
  
for i = 1:length(x)
Z(:,i) = (p_M * mean(real(V)) + q * min( max(x(i) - CM_L2, 0), CM_H2 - CM_L2 * (mean(real(V))^(2/3)) )) * c_T;
end
  figure(4)
  subplot(1, 2, 1)
  y = K2C(y);
  surf(x,y,Z)
  xlabel('CO2')
  ylabel('T')
  zlabel('pM flux')
   
  
 % acidification factor correction factor for assimilation
  CX_L = 380; % uatm,  Lower boundary for ingestion  
  CX_H = 450; % uatm,  Upper boundary for ingestion 
  c_X = (CX_H - x) / (CX_H - CX_L);
  
  p_AmT = c_T * cPar.p_Am;
  p_XmT = p_AmT / par.kap_X;
  pX = 0.4 .* p_XmT .* max(real(V)).^(2/3); 
  
  for i = 1:length(x)
Z(:,i) = par.kap_X .* pX .* c_X(i); %pA flux
  end
  figure(4)
  subplot(1, 2, 2)
  surf(x,y,Z)
  xlabel('CO2')
  ylabel('T')
  zlabel('pA flux')
    
  end
