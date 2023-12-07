% Construct the filename based on the scenario
filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';

% Load parameter values
res = load(filename);
par  = res.par;          % get primary parameters
cPar = parscomp_st(par); % get compound parameters


data.LN = [... lenght ~ Number of eggs
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
    62.04	691111.11
    ];
data.LN(:,1) = data.LN(:,1) / 10; % convert length from mm to cm
data.LN(:,1) = (data.LN(:,1) * par.del_M).^3;
data.LN(:,2) = data.LN(:,2); % Assuming two spawning events per year
units.LN   = {'cm', '#eggs'};  label.LN = {'Total length', 'Total Fecundity'};
temp.LN    = C2K(21);  units.temp.LN = 'K'; label.temp.LN = 'temperature';
bibkey.LN = 'Bilbao2010';


figure()
tbl = table(data.LN(:,1),data.LN(:,2));
mdl = fitlm(tbl,'Var2 ~ Var1');
plot(mdl)

