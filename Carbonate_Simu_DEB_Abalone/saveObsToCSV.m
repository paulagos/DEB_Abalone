%%%%%% save observables quantities form the simu structure into csc files

results  = '/Users/lagosjime/Documents/Work/UBO_Brest/Haliotis_tuberculata/Simulation_DEB_Abalone/Simulation_files/DEB_model v.1.0/Results/simu_abal_2024Mar27';

cd(results)

filenames = dir(fullfile(results, '*simu*'));

%% 1 - read data  from structure in the results folder
load(filenames(1).name)
varTypes = {'string','string', 'string','string','string','double','double','double','double','double'};
varNames = {'location','ssp','food','zoom','cohort','days','totalLength','wetWeight','fecundity','gonadWetWeight'};
nrows    = size(simu.obs,1)* size(simu.obs,3)*length(filenames);
data_DEB_Abalone = table('Size',[nrows, length(varNames)],'VariableTypes',varTypes, 'VariableNames', varNames);
row_idx = 1;

for f = 1:length(filenames)
    
    filename = filenames(f).name;
    
    load(filename) % loads simu
    
    % Get scenario parameters from simu or filename
    idx      = strfind(filename,'_');
    
    cohort   = repmat('1', length(simu.obs),1) ;
    location = filename(idx(3)+1:idx(4)-1);
    ssp      = filename(idx(2)+1:idx(3)-1);
    food     = filename(idx(1)+1:idx(2)-1);
    zoom     = filename(idx(5)+1:end-4);
    
    % get observables
    days           = simu.obs(:,1,1);
    totalLength    = simu.obs(:,3,1);
    wetWeight      = simu.obs(:,4,1);
    fecundity      = simu.obs(:,6,1);
    gonadWetWeight = simu.obs(:,7,1);
    
    numCohorts     = size(simu.obs,3);
    
    if numCohorts > 1
        for c = 2:numCohorts
            cohort         = [cohort; repmat(num2str(c), length(simu.obs),1) ];
            days           = [days; simu.obs(:,1,c)];
            totalLength    = [totalLength; simu.obs(:,3,c)];
            wetWeight      = [wetWeight; simu.obs(:,4,c)];
            fecundity      = [fecundity; simu.obs(:,6,c)];
            gonadWetWeight = [gonadWetWeight; simu.obs(:,7,c)];
            
        end
    end
    tmp = table(cohort,days,totalLength,wetWeight,fecundity,gonadWetWeight);
    
    data_DEB_Abalone(row_idx:length(days)+row_idx-1,5:end) = tmp;
    
    data_DEB_Abalone(row_idx:length(days)+row_idx-1,1) = table(repmat(location,length(days),1));
    data_DEB_Abalone(row_idx:length(days)+row_idx-1,2) = table(repmat(ssp,length(days),1));
    data_DEB_Abalone(row_idx:length(days)+row_idx-1,3) = table(repmat(food,length(days),1));
    data_DEB_Abalone(row_idx:length(days)+row_idx-1,4) = table(repmat(zoom,length(days),1));

    
    row_idx = length(days)+row_idx;
    
end
