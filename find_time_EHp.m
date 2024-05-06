function condition_EHp = find_time_EHp(E_H, E_Hp, t)
% Find the time at puberty (start of reproduction)
% (E_H >= E_Hp) from model predictions for one individual,
% with standard DEB equations (Kooijman 2010)

% 15/10/2023 - Paulo F. Lagos
%---------------------------------------------------------------

%% Capture time when E_H < E_Hp condition is met
    condition_indices = find(E_H >= E_Hp); % Find all indices where the condition is met
    if ~isempty(condition_indices)
        
        % If the condition is met at least once
        first_condition_index = condition_indices(1); % Get the index of the first met condition
        condition_EHp = t(first_condition_index); % Time in days
    else
        % If the condition is not met in this generation
        condition_EHp = 0; % Set it to NaN
        
    end
    

end
