filename = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/DEB_Haliotis_control/results_Haliotis_tuberculata.mat';
load(filename); 
[stat, txtStat] = statistics_st(metaPar.model, par); printstat_st(stat, txtStat)
