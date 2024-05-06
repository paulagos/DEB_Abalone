function exportStructToCSV(structData, fileName)
    % Extract field names
    fieldNames = fieldnames(structData);
    dataForCSV = {}; % Initialize a cell array for data
    
    % Assuming each field in structData is a scalar or a vector
    for i = 1:length(fieldNames)
        fieldName = fieldNames{i};
        fieldValue = structData.(fieldName);
        if isnumeric(fieldValue) || ischar(fieldValue)
            
            % Convert numeric arrays to strings; for other types, adjust accordingly
            fieldValueStr = num2str(fieldValue);
            dataForCSV(end+1, :) = {fieldName, fieldValueStr}; % Add field name and value as string
        end
    end
    
    % Convert the cell array to a table and write to CSV
    T = cell2table(dataForCSV, 'VariableNames', {'FieldName', 'Value'});
    writetable(T, fileName);
end