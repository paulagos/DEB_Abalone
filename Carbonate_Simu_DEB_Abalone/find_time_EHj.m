function condition_times = find_time_EHj(E_H, E_Hj, t)  
% Find the time at metamorphosis
% (E_H >= E_Hj) from model predictions for one individual,
% with standard DEB equations (Kooijman 2010)

% 12/10/2023 - Paulo F. Lagos
%---------------------------------------------------------------

%% Capture time when E_H < E_Hj condition is met
    condition_indices = find(E_H >= E_Hj); % Find all indices where the condition is met
    if ~isempty(condition_indices)
        
        % If the condition is met at least once
        first_condition_index = condition_indices(1); % Get the index of the first met condition
        condition_times = t(first_condition_index); % Time in days
    else
        % If the condition is not met in this generation
        condition_times = 0; % Set it to NaN
        
    end
    