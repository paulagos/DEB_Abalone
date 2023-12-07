function [par, metaPar, txtPar] = pars_init_Haliotis_tuberculata(metaData)
%
metaPar.model = 'abj'; 
%
%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 
%
%% core primary parameters 
par.T_A = 4176.17;    free.T_A   = 0;   units.T_A   = 'K';        label.T_A   = 'arrhenius temperature';
par.T_L = 280.15;     free.T_L   = 0;   units.T_L   = 'K';        label.T_L   = 'lower boundary tolerance range';
par.T_H = 301.15;     free.T_H   = 0;   units.T_H   = 'K';        label.T_H   = 'Upper boundary tolerance range';
par.T_AL = 13516.96;  free.T_AL  = 0;   units.T_AL  = 'K';        label.T_AL  = 'arrh. temp for lower boundary';
par.T_AH = 24570.25;  free.T_AH  = 0;   units.T_AH  = 'K';        label.T_AH  = 'arrh. temp for upper boundary';
par.z = 7.0392;       free.z     = 1;   units.z     = '-';        label.z     = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m   = 'l/d.cm^2'; label.F_m   = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.2;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.022003;     free.v     = 1;   units.v     = 'cm/d';     label.v     = 'energy conductance'; 
par.kap = 0.8981;     free.kap   = 1;   units.kap   = '-';        label.kap   = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 19.1757;    free.p_M   = 1;   units.p_M   = 'J/d.cm^3'; label.p_M   = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T   = 'J/d.cm^2'; label.p_T   = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 1;   units.k_J   = '1/d';      label.k_J   = 'maturity maint rate coefficient'; 
par.E_G = 3914.7586;  free.E_G   = 1;   units.E_G   = 'J/cm^3';   label.E_G   = '[E_G], spec cost for structure'; 
par.E_Hb = 7.838e-04; free.E_Hb  = 1;   units.E_Hb  = 'J';        label.E_Hb  = 'maturity at birth';
par.E_Hj = 9.839e-04; free.E_Hj  = 1;   units.E_Hj  = 'J';        label.E_Hj  = 'maturity at metamorphosis';
par.E_Hp = 1.896e+03; free.E_Hp  = 1;   units.E_Hp  = 'J';        label.E_Hp  = 'maturity at puberty'; 
par.h_a = 1.331e-8;   free.h_a   = 1;  units.h_a   = '1/d^2';    label.h_a   = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 1;   units.s_G   = '-';        label.s_G   = 'Gompertz stress coefficient'; 
%
%% secondary parameters 
par.del_Y = 0.3008;   free.del_Y = 0;   units.del_Y = '-';        label.del_Y = 'shape coefficient for the egg';
par.del_L = 0.502;    free.del_L = 1;   units.del_L = '-';        label.del_L = 'shape coefficient pre-metamorphosis';
par.del_M = 0.28218;  free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient post-metamorphosis'; 
par.f     = 1;        free.f     = 1;   units.f     = '-';        label.f     = 'scaled functional response for 0-var data';
par.f_ee  = 1;        free.f_ee  = 1;   units.f_ee  = '-';        label.f_ee  = 'scaled functional response for Enteropmorpha diet'; 
par.f_pp  = 1;        free.f_pp  = 1;   units.f_pp  = '-';        label.f_pp  = 'scaled functional response for Palmaria palmata diet';
par.f_sl  = 1;        free.f_sl  = 1;   units.f_sl  = '-';        label.f_sl  = 'scaled functional response for Saccorhiza polyschides diet'; 
%
%% acidification specific paramenters
par.CMH_larva = 1000; free.CMH_larva = 0;  units.CMH_larva = 'uatm';      label.CMH_larva  = 'higher pCO2 boundary for maintenace cost';          
par.CML_larva = 500;  free.CML_larva = 0;  units.CML_larva = 'uatm';      label.CML_larva  = 'lower pCO2 boundary for maintenacne cost'; 
par.CMH_adult = 1200; free.CMH_adult = 0;  units.CMH_adult = 'uatm';      label.CMH_adult  = 'higher pCO2 boundary for maintenace cost';          
par.CML_adult = 700;  free.CML_adult = 0;  units.CML_adult = 'uatm';      label.CML_adult  = 'lower pCO2 boundary for maintenacne cost'; 
par.q_larva   = 0.09; free.q_larva   = 1;  units.q_larva = 'J/cm^2 uatm^1 d^1'; label.q_larva = 'maintenace coefficient related to pCO2 for larva';
par.q_adult   = 0.09; free.q_adult   = 1;  units.q_adult = 'J/cm^2 uatm^1 d^1'; label.q_adult = 'maintenace coefficient related to pCO2 for adults';
par.kap_LR    = 0.79; free.kap_LR    = 0;  units.kap_LR  = '-';            label.kap_LR = 'yield of reproduction tissue lysis';
par.kap_LV    = 0.79; free.kap_LV    = 0;  units.kap_LV  = '-';            label.kap_LV = 'yield of structural  tissue lysis';
par.rho_V = 2.0921e+04; free.rho_V   = 0;  units.rho_V   = 'J g^-1';       label.rho_V  = 'energy content of 1g of structure'; 
%
%% shell growth parameters
% par.vA = 0;           free.vA  = 1;     units.vA = 'cm/d';        label.vA = 'coupling coefficient to assimilation';     
% par.vD = 2.37E-04;    free.vD  = 1;     units.vD = 'cm/d';        label.vD = 'coupling coefficient to dissipation';  
% par.vG = 0.01;        free.vG  = 1;     units.vG = 'cm/d';        label.vG = 'coupling coefficient to growth';
% par.d_Vs = 2.7;       free.d_Vs = 0;    units.d_Vs = 'g/cm^3';    label.d_Vs = 'specific density of calcium carbonate';
%
%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 
%
%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
