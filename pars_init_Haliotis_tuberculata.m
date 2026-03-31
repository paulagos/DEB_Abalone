function [par, metaPar, txtPar] = pars_init_Haliotis_tuberculata(metaData)

metaPar.model = 'abj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 7197;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'arrhenius temperature'; 
par.z = 3.0061;       free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 65;         free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.0256;       free.v     = 0;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.8261;     free.kap   = 0;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 17.1553;    free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 3938;       free.E_G   = 0;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 0.0002208; free.E_Hb  = 0;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 0.001507;  free.E_Hj  = 0;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metamorphosis'; 
par.E_Hp = 1.879e+03; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 4.985e-09;  free.h_a   = 0;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.01;       free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% pCO2 parameters 
par.CMH_adult = 1500;     free.CMH_adult = 0;   units.CMH_adult = 'muatm';  label.CMH_adult = 'higher pCO2 boundary for maintenace cost'; 
par.CMH_larva = 1000;     free.CMH_larva = 0;   units.CMH_larva = 'muatm';  label.CMH_larva = 'higher pCO2 boundary for maintenace cost'; 
par.CML_adult = 400;      free.CML_adult = 0;   units.CML_adult = 'muatm';  label.CML_adult = 'lower pCO2 boundary for maintenacne cost'; 
par.CML_larva = 300;      free.CML_larva = 0;   units.CML_larva = 'muatm';  label.CML_larva = 'lower pCO2 boundary for maintenacne cost'; 
par.q_adult = 0.001;      free.q_adult = 1;   units.q_adult = 'J/d.cm^2.muatm';  label.q_adult = 'maintenace coefficient related to pCO2 for adults'; 
par.q_larva = 0.001008;   free.q_larva = 0;   units.q_larva = 'J/d.cm^2.muatm';  label.q_larva = 'maintenace coefficient related to pCO2 for larva'; 

%% other parameters
par.L0_Auzous2021 = 6.023;  free.L0_Auzous2021 = 0;   units.L0_Auzous2021 = 'cm';  label.L0_Auzous2021 = 'L at start Auzous2021'; 
par.T_AH = 30068;       free.T_AH  = 0;   units.T_AH = 'K';         label.T_AH = 'arrh. temp for upper boundary'; 
par.T_AL = 908;         free.T_AL  = 0;   units.T_AL = 'K';         label.T_AL = 'arrh. temp for lower boundary'; 
par.T_H = 301.15;       free.T_H   = 0;   units.T_H = 'K';          label.T_H = 'Upper boundary tolerance range'; 
par.T_L = 279.15;       free.T_L   = 0;   units.T_L = 'K';          label.T_L = 'lower boundary tolerance range'; 
par.del_L = 0.4887;     free.del_L = 1;   units.del_L = '-';        label.del_L = 'shape coefficient pre-metamorphosis'; 
par.del_M = 0.39857;    free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient post-metamorphosis'; 
par.f = 0.9446;         free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.f_ee = 0.9;         free.f_ee  = 0;   units.f_ee = '-';         label.f_ee = 'scaled functional response for Enteropmorpha diet'; 
par.f_pp = 1;           free.f_pp  = 0;   units.f_pp = '-';         label.f_pp = 'scaled functional response for Palmaria palmata diet'; 
par.f_sl = 0.468;       free.f_sl  = 1;   units.f_sl = '-';         label.f_sl = 'scaled functional response for Saccorhiza polyschides diet'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
