function [temp, pco2] = forcing_deb(location, scenario,bdate)

% Load insitu data for all locations within france
switch location
    case 'brest'
        filename = 'Environmental_data/Brest_MAREL.csv';
    case 'roscoff'
        filename = 'Environmental_data/Roscoff_FIGUIER.csv';
    case 'nantes'
         filename = 'Environmental_data/Nantes_MOLIT.csv';
    case 'larochelle'
        filename = 'Environmental_data/Larochelle_ROCHL.csv';
    case 'med'
        filename = 'Environmental_data/Mediterranean_BESS.csv';  
end

env_data = table2array(readtable(filename));
env_data = [env_data(bdate:end,:); env_data(1:(bdate-1),:)];


% The SSP predictions for temperature and pCO2
fYears  = [2020,  2030,  2040,  2050,  2060,  2070,  2080,  2090,  2100];

ssppco2  = [417.889	451.997	491.255	535.221	583.692	636.942	695.674	760.992	834.341    % SSP3-Baseline
            412.142	424.682	425.527	420.476	413.876	406.784	398.628	389.057	377.897    % SSP1-19
            414.550	436.543	452.291	460.693	462.902	460.642	453.862	443.091	432.171    % SSP1-26
            414.723	440.921	468.530	496.066	521.324	542.585	558.132	565.688	566.859    % SSP2-45
            417.249	452.767	499.678	559.692	635.793	730.026	841.520	963.842	1088.970]; % SSP5-baseline


sspTemp = [1.226	1.524	1.864	2.203	2.551	2.912	3.284	3.672	4.071   % SSP3-Baseline
           1.237	1.466	1.558	1.562	1.528	1.493	1.456	1.410	1.345   % SSP1-19
           1.218	1.473	1.655	1.759	1.812	1.835	1.835	1.804	1.761   % SSP1-26
           1.242	1.486	1.743	1.966	2.154	2.318	2.459	2.567	2.629   % SSP2-45
           1.266	1.603	2.013	2.480	2.976	3.504	4.050	4.576	5.052]; % SSP5-baseline

% Calculate increase rates for each time interval
yearInterval = diff(fYears);

pco2Rates = diff(ssppco2, 1, 2) ./ (yearInterval*365);
tempRates = diff(sspTemp, 1, 2) ./ (yearInterval*365);

switch scenario
    case 'SSP3sBaseline'; fRow = 1;
    case 'SSP1s19';       fRow = 2;
    case 'SSP1s26';       fRow = 3;
    case 'SSP1s45';       fRow = 4;
    case 'SSP5sBaseline'; fRow = 5;
end

if strcmp(scenario, 'Current') % Get the daily temp and pco2 from 2020-2100 if running Current conditions

    temp = repmat(env_data(:,2)', 1,(fYears(end)-fYears(1)) );
    pco2 = repmat(env_data(:,3)', 1,(fYears(end)-fYears(1)) );
    
else   % Get the daily temp and pco2 from 2020-2100 if running an SSP
       
    nDays  = 365 * (fYears(end) - fYears(1));
    temp   = zeros(1,nDays);
    pco2   = temp;
    
    % Loop through each time step in the simulation
    for i = 1 : nDays
        
        % Which decade or time period are we in (i.e. which column of forcing to use)
        fCol = ceil(i / (10*365));
        periodDay = i - (fCol - 1) * 10 * 365;
        
        % Get baseline temperature and pco2 at day i
        days = 365 ;
        day = int32(((i / days) - floor(i / days)) * days);
        if (day==0); day=365; end
        T = env_data(day,2);
        P = env_data(day,3);
        
        % Get the temp and pco2 increase til day i
        cumuT = 0; cumuP = 0;
        for j = 1:fCol
            if j==fCol
                cumuT = cumuT + tempRates(fRow,j)*periodDay;
                cumuP = cumuP + pco2Rates(fRow,j)*periodDay;
            else
                cumuT = cumuT + tempRates(fRow,j)*10*365;
                cumuP = cumuP + pco2Rates(fRow,j)*10*365;
            end
        end
        
        % Get temp and pco2 at day i
        temp(i) = T + cumuT;
        pco2(i) = P + cumuP;
    end
end

temp(length(temp)+1) = temp(end); % need one more timestep bc ode starts at t=0
pco2(length(pco2)+1) = pco2(end);

end

