function [par, metaPar, txtPar] = pars_init_Haliotis_tuberculata(metaData)
%
metaPar.model = 'abj'; 
%
%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 
%
%% core primary parameters 
par.T_A = 8000;       free.T_A   = 1;   units.T_A   = 'K';        label.T_A   = 'arrhenius temperature';
par.T_L = 281;        free.T_L   = 0;   units.T_L   = 'K';        label.T_L   = 'lower boundary tolerance range';
par.T_H = 297;        free.T_H   = 0;   units.T_H   = 'K';        label.T_H   = 'Upper boundary tolerance range';
par.T_AL = 23434;     free.T_AL  = 1;   units.T_AL   = 'K';       label.T_AL  = 'arrh. temp for lower boundary';
par.T_AH = 92422;     free.T_AH  = 1;   units.T_AH  = 'K';        label.T_AH  = 'arrh. temp for upper boundary';
par.z = 3.0392;       free.z     = 1;   units.z     = '-';        label.z     = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m   = 'l/d.cm^2'; label.F_m   = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.2;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.022003;     free.v     = 1;   units.v     = 'cm/d';     label.v     = 'energy conductance'; 
par.kap = 0.86981;    free.kap   = 1;   units.kap   = '-';        label.kap   = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 19.1757;    free.p_M   = 1;   units.p_M   = 'J/d.cm^3'; label.p_M   = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T   = 'J/d.cm^2'; label.p_T   = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J   = '1/d';      label.k_J   = 'maturity maint rate coefficient'; 
par.E_G = 3914.7586;  free.E_G   = 1;   units.E_G   = 'J/cm^3';   label.E_G   = '[E_G], spec cost for structure'; 
par.E_Hb = 9.838e-04; free.E_Hb  = 1;   units.E_Hb  = 'J';        label.E_Hb  = 'maturity at birth';
par.E_Hj = 9.839e-04; free.E_Hj  = 1;   units.E_Hj  = 'J';        label.E_Hj  = 'maturity at metamorphosis';
par.E_Hp = 1.896e+03; free.E_Hp  = 1;   units.E_Hp  = 'J';        label.E_Hp  = 'maturity at puberty'; 
par.h_a = 5.331e-9;  free.h_a    = 1;   units.h_a   = '1/d^2';    label.h_a   = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G   = '-';        label.s_G   = 'Gompertz stress coefficient'; 
%
%% secondary parameters 
par.del_Y = 0.971;    free.del_Y = 1;   units.del_Y  = '-';       label.del_Y = 'shape coefficient for the egg';
par.del_L = 0.202;    free.del_L = 1;   units.del_L = '-';        label.del_L = 'shape coefficient pre-metamorphosis';
par.del_M = 0.28218;  free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient post-metamorphosis'; 
par.f     = 1;        free.f     = 1;   units.f     = '-';        label.f     = 'scaled functional response for 0-var data';
par.f_ee  = 1;        free.f_ee  = 1;   units.f_ee  = '-';        label.f_ee  = 'scaled functional response for Enteropmorpha diet'; 
par.f_pp  = 1;        free.f_pp  = 1;   units.f_pp  = '-';        label.f_pp  = 'scaled functional response for Palmaria palmata diet';
par.f_sl  = 1;        free.f_sl  = 1;   units.f_sl  = '-';        label.f_sl  = 'scaled functional response for Saccorhiza polyschides diet'; 
%
%% acidification specific paramenters
par.CX_L = 1010;      free.CX_L = 1;    units.CX_L = 'uatm';           label.CX_L   = 'lower boundary for ingestion';  
par.CX_H = 6779;      free.CX_H = 1;    units.CX_H = 'uatm';           label.CX_H   = 'upper boundary for ingestion'; 
par.CM_H = 1300;      free.CM_H = 1;    units.CM_H = 'uatm';           label.CM_H   = 'higher pCO2 boundary for maintenace cost';          
par.CM_L = 700;       free.CM_L = 1;    units.CM_L = 'uatm';           label.CM_L   = 'lower pCO2 boundary for maintenacne cost';           
par.q    = 0.009;     free.q    = 1;    units.q = 'J/cm^2 uatm^1 d^1'; label.q      = 'maintenace coefficient related to pCO2';
par.kap_LR = 0.79;    free.kap_LR = 1;  units.kap_LR = '-';            label.kap_LR = 'yield of reproduction tissue lysis';
par.kap_LV = 0.79;    free.kap_LV = 1;  units.kap_LV = '-';            label.kap_LV = 'yield of structural  tissue lysis';
par.rho_V  = 19600;   free.rho_V  = 1;  units.rho_V  = 'J g^-1';       label.rho_V  = 'energy content of 1g of structure'; 
%
%% shell growth parameters
par.vA = 0;           free.vA  = 1;     units.vA = 'cm/d';        label.vA = 'coupling coefficient to assimilation';     
par.vD = 2.37E-04;    free.vD  = 1;     units.vD = 'cm/d';        label.vD = 'coupling coefficient to dissipation';  
par.vG = 0.01;        free.vG  = 1;     units.vG = 'cm/d';        label.vG = 'coupling coefficient to growth';
par.d_Vs = 2.7;       free.d_Vs = 0;    units.d_Vs = 'g/cm^3';    label.d_Vs = 'specific density of calcium carbonate';
%
%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 
%
%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
