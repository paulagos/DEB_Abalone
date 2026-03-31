  function [data, auxData, metaData, txtData, weights] = mydata_Haliotis_tuberculata
  
%% set metaData
metaData.phylum     = 'Mollusca'; 
metaData.class      = 'Gastropoda'; 
metaData.order      = 'Vetigastropoda'; 
metaData.family     = 'Haliotidae';
metaData.species    = 'Haliotis_tuberculata'; 
metaData.species_en = 'European abalone'; 

metaData.ecoCode.climate = {'MB','MC'};
metaData.ecoCode.ecozone = {'MAE','MANE','MAm'};
metaData.ecoCode.habitat = {'0jMp', 'jiMb'};
metaData.ecoCode.embryo  = {'Mp'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'bjPp', 'jiHa'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(20); % K, body temp
metaData.data_0     = {'ab'; 'ap'; 'am'; 'Lb'; 'Lj'; 'Lp'; 'Li'; 'Ww0'; 'Wwi'}; 
metaData.data_1     = {'t-Ww'; 't-L'; 'L-Ww'; 'L-N'; 'L-WwR'}; 

metaData.COMPLETE = 2.8; % using criteria of LikaKear2011

metaData.author   = {'Paulo F. Lagos','Alva Curtsdotter','Sabine Roussel','Jonathan Flye-Sainte-Maire','Stephanie Auzoux-Bordenave','Sophie Martin','Laure Pecquerie'};    
metaData.date_subm = [2026 02 17];              
metaData.email    = {'pau.lagosj@gmail.com'};            
metaData.address  = {'LEMAR'};   

metaData.curator     = {'Bas Kooijman'};
metaData.email_cur   = {'bas.kooijman@vu.nl'}; 
metaData.date_acc    = [2026 02 17]; 

%% Zero-variate data
data.ab = 0.91;     units.ab = 'd';    label.ab = 'age at birth';               bibkey.ab = 'Wessel2018';
  temp.ab = C2K(17);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
  comment.ab = 'We assume feeding beginning at late trochophore stage (19hpf to 24hpf). However, larva might not need to eat until metamorphosis has occured (9 dpf approx)';
data.aj = 3.9;      units.aj = 'd';    label.aj = 'age at metamorphosis'; bibkey.aj = 'Auzoux2022';
  temp.aj = C2K(15);  units.temp.aj = 'K'; label.temp.aj = 'temperature';
  comment.aj = 'veliger larva at end of planktonic phase 94hpf ';
data.ap = 365*1.2;  units.ap = 'd';     label.ap = 'age since birth at puberty'; bibkey.ap = 'Coch1980';
  temp.ap = C2K(17);   units.temp.ap = 'K'; label.temp.ap = 'temperature';
  comment.ap = 'females mature at ~1 years with variations depending on locations';
data.am = 365*15; units.am = 'd';   label.am = 'life span';                bibkey.am = 'Basu1997';
  temp.am = C2K(18);  units.temp.am = 'K'; label.temp.am = 'temperature';

data.Lb = 0.0099; units.Lb = 'cm';  label.Lb = 'length at birth';          bibkey.Lb = 'Wessel2018';
  comment.Lb = 'Minimiun length for a larvae at 40hpf from experiments at 17deg C';
data.Lj = 0.0268;  units.Lj = 'cm'; label.Lj = 'length at metamorphosis';   bibkey.Lj = 'Auzoux2022';
  comment.Lj = 'Average length of 200 pre-metamorphic larva (96hph) at ph8.0 ';
data.Lp  = 3.0;  units.Lp  = 'cm';  label.Lp  = 'length at puberty';        bibkey.Lp  = 'Basu1997';
  comment.Lp = 'range from 30mm to 49mm';
data.Li  = 14; units.Li  = 'cm';  label.Li  = 'ultimate total lenght';     bibkey.Li  = 'Basu1997';
  comment.Li = 'based on one individual found on Nord-West of Cotentin (Author personal note)';

data.Ww0 = 1.77e-6;  units.Ww0 = 'g'; label.Ww0 = 'initial wet-weight';    bibkey.Ww0 = 'Bilbao2010';
  comment.Ww0 = 'based on mature oocyte size of 90-190 mum: pi/6*0.015^3';

data.Wwi = 450 * 0.73;  units.Wwi = 'g'; label.Wwi = 'maximum shell free wet-weight'; bibkey.Wwi = 'Basu1997';
  comment.Wwi = 'based on one individual found on Nord-West of Cotentin (Author personal note), 27% of the total wet-weight is the shell';

%% Uni-variate data
%
% time ~ total wet-weight
data.tW = [ ... % age (d), total wet-weight (g)
163   0.318
193   0.644
225   3.157
321   7.883
416   14.485
511   22.650
606   28.314
702   36.478
797   44.329
892   49.056
987   52.533];
data.tW(:,2) = data.tW(:,2) * 0.73; % substract 27% of the weight to account for the pressence of shell in the data
units.tW   = {'d', 'g'}; label.tW = {'Time', 'Total wet-weight'};  
temp.tW    = C2K(18); units.temp.tW = 'K'; label.temp.tW = 'temperature';
bibkey.tW  = 'Basu1997';
comment.tW = 'Age of adult abalone from the bay of Brest in months converted to days (2.7  years in total)';
%
% Data for Females: length ~ total wet-weight (including shell)
data.LWw_f = [ ... % total lenght (mm), total body wet-weight (g)
30.00	 3.56
30.56	 3.56
32.05	 4.25
33.08	 4.71
34.01	 4.94
35.04	 5.17
35.04	 5.86
36.06	 7.47
36.44	 5.40
37.00	 6.78
36.90	 6.09
37.56	 5.86
37.56	 7.01
37.93	 6.55
38.58	 7.24
38.58	 6.32
38.96	 7.24
39.42	 8.16
40.07	 8.39
39.98	 7.24
40.45	 7.24
40.45	 8.16
40.82	 9.08
40.91	 8.39
40.91	 7.70
41.47	 8.39
42.03	 7.24
41.85	 8.62
41.94	 10.23
41.94	 10.92
42.59	 9.31
43.06	 11.61
42.97	 10.23
42.97	 8.62
43.34	 10.46
44.09	 11.61
43.99	 10.69
44.09	 9.77
44.46	 10.00
44.55	 11.15
45.02	 13.22
45.02	 12.53
44.93	 11.38
44.93	 10.46
45.02	 9.54
45.58	 10.69
45.58	 13.45
45.86	 13.68
46.04	 12.53
46.04	 11.38
46.04	 10.46
46.04	 15.98
46.98	 14.14
46.98	 13.45
47.07	 12.53
47.07	 11.61
47.54	 10.46
48.94	 7.47
49.50	 10.46
48.10	 11.84
48.10	 12.99
48.19	 14.83
48.19	 15.98
48.00	 16.67
48.56	 17.36
47.44	 15.98
47.44	 15.52
48.94	 16.44
48.38	 14.37
48.94	 14.37
49.50	 14.14
50.43	 13.68
50.43	 15.75
50.06	 15.29
49.50	 15.29
49.68	 16.44
48.84	 15.52
49.96	 17.13
49.96	 17.59
49.96	 18.97
49.96	 20.11
51.08	 23.33
51.46	 20.80
51.08	 19.66
50.99	 18.97
51.08	 17.13
51.92	 14.83
51.92	 15.52
52.01	 17.13
51.92	 18.05
52.95	 18.51
53.04	 19.66
52.95	 20.80
54.07	 21.72
54.07	 24.02
53.97	 23.79
55.09	 24.02
55.00	 22.18
54.91	 19.66
55.00	 18.97
55.93	 19.66
55.93	 20.80
55.93	 22.41
55.47	 21.95
56.03	 24.25
55.93	 25.86
55.93	 27.24
55.93	 30.69
56.96	 30.69
57.52	 24.71
59.01	 23.79
59.48	 28.39
61.06	 31.84
63.12	 30.92
64.61	 32.30
63.02	 34.60
60.97	 38.28
64.98	 48.16
64.98	 39.66
66.01	 39.66
71.98	 52.99
71.42	 58.05];
data.LWw_f(:,1) = sort(data.LWw_f(:,1)/10); %convert length from mm to cm
data.LWw_f(:,2) = data.LWw_f(:,2) * 0.73; % the weight of the shell is 27% of total wet-weight
units.LWw_f   = {'cm', 'g'};  label.LWw_f = {'Total length', 'Total wet-weight'};  
%temp.LWw_f    = C2K(20);  units.temp.LWw_f = 'K'; label.temp.LWw_f = 'temperature';
bibkey.LWw_f = 'Bilbao2010';
comment.LWw_f = 'Data of abalone females from the canary islands, total body weight including shell';
%
% Length (mm) ~ fecundity (#eggs)
data.LN = [...
36.04	53333.33
38.95	87037.04
40.97	151851.85
41.02	115555.56
42.99	188148.15
44.97	182962.96
45.95	245185.19
47.98	177777.78
52.12	328148.15
56.02	312592.59
57.01	507037.04
59.03	595185.19
59.08	634074.07
62.04	691111.11];
data.LN(:,1) = data.LN(:,1)/10; % Convert length from mm to cm
units.LN   = {'cm', '#'};  label.LN = {'Lenght', 'Fecundity'};  
temp.LN    = C2K(20);  units.temp.LN = 'K'; label.temp.LN = 'temperature';
bibkey.LN  = 'Bilbao2010';
comment.LN = 'Batch fecundity for animals from the Canary islands estimated by counting the number of stage III and IV eggs in wild females';
%
% Total lenght (mm) - Total weight (g)
data_LWw_Roussel2019 = readmatrix('data_LWw_Roussel2019.txt');
data_LWw_Roussel2019(:,1) =  data_LWw_Roussel2019(:,1) / 10;  % convert length from mm to cm
data_LWw_Roussel2019(:,2) = data_LWw_Roussel2019(:,2) * 0.73;  % loss 1/3 of the weight to account for the pressence of shell in the data
[L,~,ci] = unique(data_LWw_Roussel2019(:,1));
maxci     = max(ci);
Ww_mean     = zeros(size(L,1),1);
Ww_sd       = zeros(size(L,1),1);
for i = 1:maxci
    Ww_mean(i) = mean(data_LWw_Roussel2019(ci == i,2));
    Ww_sd(i)   = std(data_LWw_Roussel2019(ci == i,2)); 
end
data.LWw_Roussel2019 = [L, Ww_mean]; %data  Total lenght (cm) - Total weight (g)
units.LWw_Roussel2019 = {'cm', 'g'};  label.LWw_Roussel2019 = {'Total length', 'Total weight'};  
%temp.LWw_Roussel2019 = C2K(21);  units.temp.LWw_Roussel2019 = 'K'; label.temp.LWw_Roussel2019 = 'temperature';
stdev.LWw_Roussel2019 = Ww_sd ; units.stdev.LWw_Roussel2019 = 'cm'; label.stdev.LWw_Roussel2019 = 'standard deviation';
bibkey.LWw_Roussel2019 = 'Roussel2019';
comment.LWw_Roussel2019 = '10 to 34 months old farmed abalone from france Haliotis';
%
% Roussel (2019)
% Data for females under a diet of Saccorhiza polyschides
M = 1; F = 2;
Roussel2019_sl = [... % (1)total lenght(cm)~ (2)sex ~ (3)gonad weight(g) ~ (4) total weight(g)
3.45	F	0.68	7.93
3.20	M	0.52	6.42
3.88	F	0.58	10.24
4.34	F	1.37	16.47
3.12	F	0.28	5.79
4.50	M	1.37	18.61
3.17	F	0.37	4.84
2.92	F	0.29	4.47
2.57	F	0.26	3.3
3.07	F	0.17	4.66
3.20	M	0.37	6.16
3.87	F	0.42	10.21
3.66	F	0.64	9.53
3.55	F	0.39	7.76
3.30	F	0.32	5.82
3.92	M	0.90	11.75
3.52	F	0.32	7.13
3.22	F	0.38	5.95
3.10	F	0.31	5.28
2.97	F	0.35	4.64
3.60	M	0.64	8.10
2.65	M	0.24	3.48
3.10	M	0.42	5.24
3.33	F	0.54	7.71
3.00	F	0.45	5.81
3.53	M	0.67	8.90
3.83	M	0.67	10.37
3.35	M	0.39	6.86
3.20	M	0.28	5.76
3.40	M	0.46	7.61
3.00	M	0.31	5.01];
% Subtract 1/3 of the total weight (4th column) to account for the shell
Roussel2019_sl(:,4) = Roussel2019_sl(:,4) * 0.73;

% Separate female data
female_data_sl = Roussel2019_sl(Roussel2019_sl(:,2) == 2, :);
%
data.LGw_sl_female = sort(female_data_sl(:,[1 3])); % Female total length(cm) ~ gonad weight(g)
units.LGw_sl_female = {'cm', 'g'}; label.LGw_sl_female = {'total length', 'gonad weight'};  
temp.LGw_sl_female = C2K(14.4); units.temp.LGw_sl_female = 'K'; label.temp.LGw_sl_female = 'temperature';
bibkey.LGw_sl_female = 'Roussel2019';
comment.LGw_sl_female = 'one year old adult female abalone fed for 11 months with a diet of Saccorhiza polyschides ';
%
data.LWw_sl_female = sort(female_data_sl(:,[1 4])); % Female total length(cm) ~ total weight(g)
units.LWw_sl_female = {'cm', 'g'}; label.LWw_sl_female = {'total length', 'wet-weight'};  
%temp.LWw_sl_female = C2K(14.4); units.temp.LWw_sl_female = 'K'; label.temp.LWw_sl_female = 'temperature';
bibkey.LWw_sl_female = 'Roussel2019';
comment.LWw_sl_female = 'one year old adult female abalone fed for 11 months with a diet of Saccorhiza polyschides';
%
% Data for male and female under a diet of Palmaria palmata
Roussel2019_pp = [... % (1)total length(cm)~ (2)sex ~ (3)gonad weight(g) ~ (4) total weight(g)
4.30	F	1.02	15.88
3.37	M	0.67	8.270
3.23	F	0.48	8.320
3.88	F	1.07	14.33
3.60	F	0.67	11.25
3.75	F	0.79	10.74
4.24	F	1.92	21.31
3.68	F	0.78	10.95
4.16	F	1.37	17.03
3.90	F	0.82	13.27
3.86	M	1.05	13.40
3.30	M	0.62	7.470
3.63	F	0.96	9.920
3.50	F	0.55	8.560
4.25	F	1.32	15.40
3.15	M	0.49	6.740
3.60	F	0.84	10.51
3.83	F	0.71	12.56
3.31	M	0.57	8.230
3.95	F	1.00	14.66
4.15	M	1.50	14.24
4.12	M	1.28	15.06
3.55	M	0.70	9.380
3.67	M	0.72	9.650
4.17	F	1.20	18.49
3.64	M	0.78	10.28
3.56	M	0.54	9.660
4.20	M	1.13	14.38
5.15	F	2.17	31.47
3.77	F	0.61	10.99
4.60	F	1.82	21.18
4.70	F	2.01	23.42
4.10	F	1.06	16.78
4.10	F	1.33	15.24
4.70	M	1.76	19.35
3.83	M	0.83	13.00];
% Subtract 1/3 of the total weight (4th column) to account for the shell
Roussel2019_pp(:,4) = Roussel2019_pp(:,4) * 0.73;

% Separate data by sex
female_data_pp = Roussel2019_pp(Roussel2019_pp(:,2) == 2, :);
%
data.LGw_pp_female = sort(female_data_pp(:,[1 3])); % total length(cm) ~ gonad weight(g)
units.LGw_pp_female = {'cm', 'g'}; label.LGw_pp_female = {'total length', 'gonad weight'};  
temp.LGw_pp_female = C2K(14.4); units.temp.LGw_pp_female = 'K'; label.temp.LGw_pp_female = 'temperature';
bibkey.LGw_pp_female = 'Roussel2019';
comment.LGw_pp_female = 'one year old juveniles female abalone fed for 11 months with a diet of Palmaria palmata';
%
data.LWw_pp_female = sort(female_data_pp(:,[1 4])); % total length(cm) ~ total weight(g)
units.LWw_pp_female = {'cm', 'g'}; label.LWw_pp_female = {'total length', 'wet-weight'};  
temp.LWw_pp_female = C2K(14.4); units.temp.LWw_pp_female = 'K'; label.temp.LWw_pp_female = 'temperature';
bibkey.LWw_pp_female = 'Roussel2019';
comment.LWw_pp_female = 'one year old juveniles female abalone fed for 11 months with a diet of palmaria palmata';
%
% Data for male and female under a diet of Enteromorpha sp.
Roussel2019_ee = [... % (1)total length(cm)~ (2)sex ~ (3)gonad weight(g) ~ (4) total weight(g)
3.56	F	1.04	9.04
3.63	M	0.66	8.25
3.52	M	1.01	9.14
4.05	F	1.05	11.71
3.30	F	0.68	7.64
3.40	M	0.63	7.04
3.63	F	0.73	8.55
3.34	M	0.70	6.63
3.40	M	0.56	7.34
3.85	F	0.93	10.22
3.60	F	0.92	9.74
3.25	F	0.79	7.06
2.94	M	0.48	4.75
3.57	F	0.81	9.3
3.27	F	0.68	7.23
3.18	F	0.44	5.93
3.50	F	0.70	7.88
3.24	M	0.62	5.95
3.70	F	0.65	8.59
3.68	M	0.76	9.13
4.15	M	0.65	12.43
3.35	M	0.54	6.35
3.26	M	0.42	6.31
4.05	M	0.90	11.95
3.48	F	0.47	7.35
3.70	F	0.73	8.64
3.13	M	0.62	4.9
3.20	F	0.43	6.07
3.20	M	0.38	5.22
3.20	M	0.35	5.46
2.84	F	0.29	4.02
2.86	M	0.30	3.93
2.95	M	0.39	4.42];
% Subtract 1/3 of the total weight (4th column) to account for the shell
Roussel2019_ee(:,4) = Roussel2019_ee(:,4) * 0.73;

% Separate data by sex (1 = male, 2 = female)
female_data_ee = Roussel2019_ee(Roussel2019_ee(:,2) == 2, :);
%
data.LGw_ee_female = sort(female_data_ee(:,[1 3])); % total length(cm) ~ gonad weight(g)
units.LGw_ee_female = {'cm', 'g'}; label.LGw_ee_female = {'total length', 'gonad weight'};  
temp.LGw_ee_female = C2K(14.4); units.temp.LGw_ee_female = 'K'; label.temp.LGw_ee_female = 'temperature';
bibkey.LGw_ee_female = 'Roussel2019';
comment.LGw_ee_female = 'one year old juveniles female abalone fed for 11 months with a diet of Enteromorpha sp.';
%
data.LWw_ee_female = sort(female_data_ee(:,[1 4])); % total length(cm) ~ total weight(g)
units.LWw_ee_female = {'cm', 'g'}; label.LWw_ee_female = {'total length', 'wet-weight'};  
temp.LWw_ee_female = C2K(14.4); units.temp.LWw_ee_female = 'K'; label.temp.LWw_ee = 'temperature';
bibkey.LWw_ee_female = 'Roussel2019';
comment.LWw_ee_female = 'one year old juveniles female abalone fed for 11 months with a diet of Enteromorpha sp.';
%
% temperatures at 3mt depth for 15 days starting the 9th of November 2021. Obtained from the COAST-HF-Iroise buoy time series for the Bay of Brest
tT = [...% time (d), temperature (°C)
0	14.04
1	13.96
2	13.95
3	13.95
4	13.88
5	13.82
6	13.80
7	13.79
8	13.72
9	13.66
10	13.54
11	13.48
12	13.25
13	13.01
14	13.87];
tT(:,1) = tT(:,1) + 1096; % convert days to match experimental time 
tT = [0, 14.04; ...
      tT; ...
      15 + 1096, 13.87];
tT(:,2) = C2K(tT(:,2)); % convert temperatute from (°C) to (°K)
%
% time (d) ~ total length (mm)
data_tL_Auzoux2021 = readmatrix("data_tL_Auzoux2021.txt");
data_tL_Auzoux2021(:,2) = data_tL_Auzoux2021(:,2) / 10; % convert mm to cm
[t,~,ci] = unique(data_tL_Auzoux2021(:,1));
maxci     = max(ci);
Lmean     = zeros(size(t,1),1);
Lsd       = zeros(size(t,1),1);
for i = 1:maxci
    Lmean(i) = mean(data_tL_Auzoux2021 (ci == i,2));
    Lsd(i)   = std(data_tL_Auzoux2021 (ci == i,2)); 
end
data.tL_auzoux2021 = [t, Lmean];
units.tL_auzoux2021 = {'d', 'cm'}; label.tL_auzoux2021 = {'time', 'total length'};  
temp.tL_auzoux2021 = tT; units.temp.tL_auzoux2021 = {'d','K'}; label.temp.tL_auzoux2021 = {'days','temperature'};
pCO2.tL_auzoux2021 = 567; units.pCO2.tL_auzoux2021 = 'μatm'; label.pCO2.tL_auzoux2021 = 'pCO2';
stdev.tL_auzoux2021 = Lsd ; units.stdev.tL_auzoux2021 = 'cm'; label.stdev.tL_auzoux2021 = 'standard deviation';
bibkey.tL_auzoux2021 = 'Auzoux2021';
comment.tL_auzoux2021 = 'Control treatment ambient pCO2 = 567μatm (pH = 7.9) three year old adult abalone fed ad libidum for 15 days from France Haliotis sea cages';

% Control: time (d) ~ larva length (mm)
data_tL_Auzoux2022 = readmatrix("data_tL_azoux2022.txt");
data_tL_Auzoux2022(:,1) = data_tL_Auzoux2022(:,1) / 24; % convert hpf to days
data_tL_Auzoux2022(:,2) = data_tL_Auzoux2022(:,2) / 10; % convert mm to cm
[t,~,ci] = unique(data_tL_Auzoux2022(:,1));
maxci     = max(ci);
Lmean     = zeros(size(t,1),1);
Lsd       = zeros(size(t,1),1);
for i = 1:maxci
    Lmean(i) = mean(data_tL_Auzoux2022 (ci == i,2));
    Lsd(i)   = std(data_tL_Auzoux2022 (ci == i,2)); 
end
data.tL_auzoux2022 = [t, Lmean];
units.tL_auzoux2022 = {'d', 'cm'}; label.tL_auzoux2022 = {'time', 'Larval length'};  
temp.tL_auzoux2022 = C2K(19.2); units.temp.tL_auzoux2022 = {'d','K'}; label.temp.tL_auzoux2022 = {'days','temperature'};
pCO2.tL_auzoux2022 = 439; units.pCO2.tL_auzoux2022 = 'μatm'; label.pCO2.tL_auzoux2022 = 'pCO2';
stdev.tL_auzoux2022 = Lsd ; units.stdev.tL_auzoux2022 = 'cm'; label.stdev.tL_auzoux2022 = 'standard deviation';
bibkey.tL_auzoux2022 = 'Auzoux2022';
comment.tL_auzoux2022 = 'Control treatment for larval growth experiments at pCO2 = 439 μatm (pH = 8.01)';
%    
% Experimental: time (d) ~ larva length (mm)
data_tL_Auzoux2022_exp = readmatrix("PCO2_Auzoux2022.txt");
data_tL_Auzoux2022_exp(:,1) = data_tL_Auzoux2022_exp(:,1) / 24; % convert hpf to days
data_tL_Auzoux2022_exp(:,2) = data_tL_Auzoux2022_exp(:,2) / 10; % convert mm to cm
[t,~,ci] = unique(data_tL_Auzoux2022_exp(:,1));
maxci     = max(ci);
Lmean     = zeros(size(t,1),1);
Lsd       = zeros(size(t,1),1);
for i = 1:maxci
    Lmean(i) = mean(data_tL_Auzoux2022_exp (ci == i,2));
    Lsd(i)   = std(data_tL_Auzoux2022_exp (ci == i,2)); 
end
data.tL_auzoux2022_pCO2 = [t, Lmean];
units.tL_auzoux2022_pCO2 = {'d', 'cm'}; label.tL_auzoux2022_pCO2 = {'time', 'Larval length'};  
temp.tL_auzoux2022_pCO2 = C2K(19.2); units.temp.tL_auzoux2022_pCO2 = {'d','K'}; label.temp.tL_auzoux2022_pCO2 = {'days','temperature'};
pCO2.tL_auzoux2022_pCO2 = 951; units.pCO2.tL_auzoux2022_pCO2 = 'μatm'; label.pCO2.tL_auzoux2022_pCO2 = 'pCO2';
stdev.tL_auzoux2022_pCO2 = Lsd ; units.stdev.tL_auzoux2022_pCO2 = 'cm'; label.stdev.tL_auzoux2022_pCO2 = 'standard deviation';
bibkey.tL_auzoux2022_pCO2 = 'Auzoux2022';
comment.tL_auzoux2022_pCO2 = 'High PCO2 treatment for larval growth experiments';

% Shell lenght (cm) ~ total wet-weight (g) 
data_LW_Chavaud2019 = readmatrix('data_LWw_Chavaud2019.txt');
data_LW_Chavaud2019(:,1) = sort(data_LW_Chavaud2019(:,1));
data_LW_Chavaud2019(:,2) = data_LW_Chavaud2019(:,2) * 0.73; % loss 1/3 of the weight to account for the pressence of shell in the data
data.LW_chauvaud = data_LW_Chavaud2019;
units.LW_chauvaud = {'cm', 'g'}; label.LW_chauvaud = {'total length', 'total weight'};  
temp.LW_chauvaud = C2K(11.5); units.temp.LW_chauvaud = {'K'}; label.temp.LW_chauvaud = {'temperature'};
pCO2.LW_chauvaud = 400; units.pCO2.LW_chauvaud = 'μatm'; label.pCO2.LW_chauvaud = 'pCO2';
bibkey.LW_chauvaud = 'Chauvaud2022';
comment.LW_chauvaud = 'Wild abalone from the bay of Brest from march 2018';

%% set weights for all real data
weights = setweights(data, []);

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
weights.psd.k_J = 0; weights.psd.k = 0.1;

%% pack auxData and txtData for output
auxData.temp    = temp;
auxData.pCO2    = pCO2;
auxData.stdev   = stdev;
txtData.units   = units;
txtData.label   = label;
txtData.bibkey  = bibkey;
txtData.comment = comment;
%
%% Group plots 
set1 = {'LGw_ee_female','LGw_pp_female','LGw_sl_female'}; subtitle1 = {'Gonad weight at several diets for female abalone'};
set2 = {'LWw_sl_female','LWw_pp_female','LWw_ee_female'}; subtitle2 = {'Lenght~weigth at several diets for female abalone'};
set3 ={'tL_auzoux2022','tL_auzoux2022_pCO2'}; subtitle3 = {'Control and experimantal pCO2 data'};
metaData.grp.sets = {set1, set2, set3};            
metaData.grp.subtitle = {subtitle1,subtitle2,subtitle3};

%% Discussion points
D1 = 'We assumed that the weight of the shell is 25% of the total wet-weight of the animals, based on calculations done on 150 abalone of aprox 2 years of age';
D2 = 'Larval stages can eat with feeding starting after metamorphosis';
metaData.discussion = struct('D2', D2);

%% Links
metaData.links.id_CoL   = '6KZCH'; % Catalogue of Life
metaData.links.id_ITIS  = '69508'; % ITIS
metaData.links.id_EoL   = '620404'; % Encyclopedia of Life
metaData.links.id_Wiki  = 'Haliotis_tuberculata'; % Wikipedia
metaData.links.id_ADW   = 'Haliotis_tuberculata'; % Animal deiversity web
metaData.links.id_WoRMS = '140059'; % WoRMS
metaData.links.id_molluscabase = '140059'; % molluscabase

%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{https://en.wikipedia.org/wiki/Haliotis_tuberculata}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{../../../bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Coch1980'; type = 'phdthesis'; bib = [ ...
'author = {Cochard, J}, ' ...
'year = {1980}, ' ...
'title  = {Recherches sur les facteurs determinant al sexualit\''{e} et la reproduction chez \emph{Haliotis tuberculata} {L}}, ' ...
'school = {Universit\''{e} de Bretagne Occidentale}, ' ...
'pages = {1-163},'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Auzoux2022'; type = 'Article'; bib = [ ... 
'author = {﻿Auzoux-Bordenave, S. and Ledoux, A. and Martin, S. and Di Poi, C. and  Suquet, M. and Badou, A. and Gaillard, F. and Servili, A. and Le Goïc, N. and Huchette, S. and Roussel, S.},'...
'year = {2022},' ...
'title = {Responses of early life stages of {E}uropean abalone (\emph{Haliotis tuberculata}) to ocean acidification after parental conditioning: {I}nsights from a transgenerational experiment}, '...
'journal = {Marine Environmental Research}, ' ...
'volume = {181}, ' ...
'pages = {1-14}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Basu1997'; type = 'Article'; bib = [ ... 
'author = {Olivier Basuyaux}, ' ... 
'year = {1997}, ' ...
'title = {Etude et mod\''{e}lisation des param\`{e}tres physico-chemiques sur la croissance de l''ormeau (\emph{Haliotis tuberculata}) en \''{e}levage en circuit semilogx-ferm\''{e}}, '...
'journal = {Fisheries Research}, ' ...
'volume = {112}, ' ...
'pages = {13-21}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Bilbao2010'; type = 'Article'; bib = [... 
'author = {Haywood, G.J. and Burns, C.W.},'...
'year = {2010},'...
'title = {Reproduction, Fecundity, and Growth of Abalone \emph{Haliotis Tuberculata Coccinea}, {R}eeve 1846 in the {C}anary {I}slands},'...
'journal = {Journal of Shellfish Research},'...
'pages = {959-967},'...
'volume = {29}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Roussel2019'; type = 'Aricle'; bib = [... 
'author = {Roussel, S. and Bisch, T. and Lachambre, S. and Boudry, P. and Gervois, J. and Lambert, C. and Huchette, S. and Day, R.},'...
'year = {2019},'...
'title = {Anti-predator response of \emph{Haliotis tuberculata} is modified after only one generation of domestication},'...
'journal = {Aquaculture Environment Interactions},'...
'volume = {11},'...
'pages = {129-142}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Wessel2018'; type = 'Article'; bib = [... 
'author = {Wessel, N. and Martin, S. and Badou, A. and Dubois, P. and Huchette, S. and Julia, V. and Nunes, F. and Harney, E. and Paillard, C. and Auzoux-Bordenave, S.},'...
'year = {2018},'...
'title = {Effect of {C}{O}2–induced ocean acidification on the early development and shell mineralization of the {E}uropean abalone (\emph{Haliotis tuberculata})},'...
'journal = {Journal of Experimental Marine Biology and Ecology},'...
'volume = {508},'...
'pages = {52-63}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Chauvaud2022'; type = 'phdthesis'; bib = [... 
'author = {Chauvaud, P.},'...
'year = {2022},'...  
'title = {Semer l'' ormeau dans le milieu naturel : r\^{o}les du milieu et des caract\''{e}ristiques biologiques d'' \emph{Haliotis tuberculata} dans la r\''{e}ussite du repeuplement},'...
'school = {Laboratoire des Sciences de l'' Environnement Marin (LEMAR)}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Auzoux2021'; type = 'Article'; bib = [... 
'author = {Auzoux-Bordenave, S. and Chevret, S. and Badou., A. and Martin, S. and Di Giglio, S. and Dubois, P.},'...
'year = {2021},'...  
'title = {Acid–base balance in the h{\ae}molymph of {E}uropean abalone (\emph{Haliotis tuberculata}) exposed to {C}{O}2-induced ocean acidification},'...
'journal = {Comparative Biochemestry and Physiology, Part A},'...
'volume = {259},'...
'pages = {1-9}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'sealifebase'; type = 'Misc'; bib = ...
'howpublished = {\url{https://www.sealifebase.ca/summary/Haliotis-tuberculata.html}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

