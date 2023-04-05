function get_plots(simu)

  par  = simu.par; % parameters
  E_Hp = par(14);
  K    = par(22); 
  Em   = par(5) / par(6); % J/cm^3  [Em] = p_Am / v

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
  
  f = X ./ (X + K); %scaled functional response
  
  % state variables
  E   = simu.tEVHRSO(:,2); % reserves
  V   = simu.tEVHRSO(:,3); % structural volume
  E_H = simu.tEVHRSO(:,4); % maturity
  E_R = simu.tEVHRSO(:,5); % reproduction buffer
  E_S = simu.tEVHRSO(:,6); % shell volume
  O   = simu.tEVHRSO(:,7); % shell opacity
  
  sc_res_dens   = E ./ (V * Em); % -, scaled reserve density e
  
  % spawning date indices
  i_sp = find(and((E_R == 0),(E_H >= E_Hp)));
  i_sp = i_sp - 1 ; % look at the preceding line with E_R value before spawning
  
   
  % observable quantities
  Lw  = simu.obs(:,2);
  Ww  = simu.obs(:,3);
  Ew  = simu.obs(:,4);
  F   = simu.obs(:,5);
  S_l = simu.obs(:,6); % shelll length
  
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
  