classdef AbaloneApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                    matlab.ui.Figure
        RunSimulationButton         matlab.ui.control.Button
        InstructionsButton          matlab.ui.control.Button
        ScenarioListBox             matlab.ui.control.ListBox
        LocationListBox             matlab.ui.control.ListBox
        FoodQualityListBox          matlab.ui.control.ListBox
        EnvVariabilityListBox       matlab.ui.control.ListBox
        ZoomFactorListBox           matlab.ui.control.TextArea
        SelectedScenariosLabel      matlab.ui.control.TextArea
        SelectedLocationsLabel      matlab.ui.control.TextArea
        SelectedFoodQualityLabel    matlab.ui.control.TextArea
        SelectedEnvVariabilityLabel matlab.ui.control.TextArea
        SelectedZoomFactorLabel     matlab.ui.control.TextArea
        GenNumberNumericField       matlab.ui.control.NumericEditField
        BornDateNumericField        matlab.ui.control.NumericEditField
    end

    methods (Access = private)
        
        % Function to update selected scenarios display
        function updateSelectedScenarios(app, event)
            selectedScenarios = app.ScenarioListBox.Value;
            app.SelectedScenariosLabel.Value = strjoin(selectedScenarios, ', ');
        end

        % Function to update selected locations display
        function updateSelectedLocations(app, event)
            selectedLocations = app.LocationListBox.Value;
            app.SelectedLocationsLabel.Value = strjoin(selectedLocations, ', ');
        end

        % Function to update selected food quality display
        function updateSelectedFoodQuality(app, event)
            selectedFoodQuality = app.FoodQualityListBox.Value;
            app.SelectedFoodQualityLabel.Value = strjoin(selectedFoodQuality, ', ');
        end

        % Function to update selected environmental variability display
        function updateSelectedEnvVariability(app, event)
            selectedEnvVariability = app.EnvVariabilityListBox.Value;
            app.SelectedEnvVariabilityLabel.Value = strjoin(selectedEnvVariability, ', ');
        end
        
        % Function to update selected zoom factor display
        function updateSelectedZoomFactor(app, event)
            zoomFactorInput = app.ZoomFactorListBox.Value;
            
            % Split the input string by commas and convert to numeric values
            zoomFactors = str2double(strsplit(zoomFactorInput{1}, ','));
            
            % Check for any invalid (NaN) inputs
            if any(isnan(zoomFactors))
                % Display an alert if there are invalid numeric values
                uialert(app.UIFigure, 'Please enter valid numeric values separated by commas.', 'Input Error');
            else
                % Convert numeric array to string array and update the label
                zoomFactorsStr = arrayfun(@num2str, zoomFactors, 'UniformOutput', false);
                app.SelectedZoomFactorLabel.Value = strjoin(zoomFactorsStr, ', ');
            end
        end

        % Callback for Run Simulation Button
        function onRunSimulationButtonPressed(app, event)
            try
                % Disable button to prevent multiple clicks
                app.RunSimulationButton.Enable = 'off';
                
                % Change pointer to hourglass (watch)
                app.UIFigure.Pointer = 'watch';
                
                % Collect user input values from UI components
                scenarios = app.ScenarioListBox.Value;
                locations = app.LocationListBox.Value;
                foodQuality = app.FoodQualityListBox.Value;
                envVariability = app.EnvVariabilityListBox.Value;
                zoomFactor = app.ZoomFactorListBox.Value;
                genNumber = app.GenNumberNumericField.Value;
                bornDate = app.BornDateNumericField.Value;
                
                % Parse zoomFactor to a numeric array
                zoomFactorParsed = str2double(strsplit(zoomFactor{1}, ','));
                
                % Assign the simulation options to the base workspace
                assignin('base', 'selectedScenarios', scenarios);
                assignin('base', 'selectedLocations', locations);
                assignin('base', 'selectedFoodQuality', foodQuality);
                assignin('base', 'selectedEnvVariability', envVariability);
                assignin('base', 'zoomFactor', zoomFactorParsed);
                assignin('base', 'genNumber', genNumber);
                assignin('base', 'bornDate', bornDate);

                % Call main.m as a script (no inputs or outputs)
                run('main.m');  % Make sure main.m is accessible in the path
                
                % Optionally, display a success message
                uialert(app.UIFigure, 'Simulation completed successfully.', 'Success');
                
                % Revert pointer to arrow after simulation
                app.UIFigure.Pointer = 'arrow';
                
                % Enable button after simulation
                app.RunSimulationButton.Enable = 'on';
                
            catch ME
                % Handle errors gracefully
                app.UIFigure.Pointer = 'arrow'; % Reset pointer in case of error
                uialert(app.UIFigure, ['Error: ', ME.message], 'Simulation Error');
                app.RunSimulationButton.Enable = 'on';
            end
        end

        % Callback for Instructions Button
        function onInstructionsButtonPressed(app, event)
            instructionFig = uifigure('Name', 'Instructions', 'Position', [400 400 500 600]);
            uilabel(instructionFig, 'Text', 'Instructions for using the Abalone DEB App:', ...
                'Position', [20 550 310 50], 'FontSize', 14, 'FontWeight', 'bold');
            
            % Add instruction labels
            uilabel(instructionFig, 'Text', {'- 1) Select the simulation scenarios.', ...
                'The simulation included five SSPs scenarios that modify key environmental',...
                'variables, such as temperature and pCO2. You can select one or multple scenarios ',...
                'from the list to pruduce several predictions simulataneusly'}, ...
                'Position', [20 495 450 70]);
            uilabel(instructionFig, 'Text', {'- 2) Choose the locations from the list.', ...
                'Currently, the simulation produce outups for five coastal locationsalong along',...
                'the west entrance to the English channel and the West Atlantic coast of France',...
                'and one location od the meditterranean sea. Using temperature and pCO2',... 
                'time series obtained from the monitoring dataset (cocoriCO2)'},...
                'Position', [20 420 450 70]);
            uilabel(instructionFig, 'Text', {'- 3) Input food quality and environmental variability.', ...
                'We use six year time series of algal boiomass from the bay of Kerazan located in',...
                'the north of Finistère, France to calculate a porcentual decline on food quality.',...
                'selected values correspond to a deacrease in food quality of 15% and 30% ',...
                'compared to normal food quality leveles 0%'},...
                'Position', [20 340 450 70]);
            uilabel(instructionFig, 'Text', {' -4) Enter number of individuales to simulate.', ...
                'The number of individuals can range from 0.5 to a maximum of 2.',...
                'To ensure the correct order of the simulation outputs, always',...
                'Input a value of 1 first, followed by other choosen values. Each value',...
                'value represent one individual'},...
                'Position', [20 260 450 70]);
            uilabel(instructionFig, 'Text', {'- 5) Number of simulations', ...
                'Simulations can be extended indefinitely. Running five simualtions extends',...
                'preditions to the end of the century. Each extention represent a period',...
                'of 20 years, corresponding to the maximun life span of H.tuberculata.'},...
                'Position', [20 190 450 70]);
            uilabel(instructionFig, 'Text', {'- 6) Start day', ...
                'The start day is currently linked to the environmental datasets in use',...
                'We do not recomend changind the start day unless new dataset are',...
                'provided. For further details contact the any of the authors.'},...
                'Position', [20 125 450 70]);
              uilabel(instructionFig, 'Text', '- 7) Run the simulation!!!', ...
                'Position', [20 80 450 70]);
            uibutton(instructionFig, 'Text', 'OK', 'Position', [150 30 100 30], ...
                'ButtonPushedFcn', @(~, ~) close(instructionFig)); % Close the window when OK is pressed
        end
    end

    % Create UI components
    methods (Access = private)

        % Code to create app components
        function createComponents(app)
            % Create UIFigure and components
            app.UIFigure = uifigure('Position', [100 100 500 700]);
            
            % Add Left EU Logo (top-left corner)
            uiimage(app.UIFigure, 'ImageSource', 'EU logo.jpg', 'Position', [10 650 100 50]);
            % Add Right Region Bretagne Logo (top-right corner)
            uiimage(app.UIFigure, 'ImageSource', 'Région Bretagne logo.jpg', 'Position', [390 650 100 50]);
            
            % Add Central Right UBO Logo (top-right corner)
            uiimage(app.UIFigure, 'ImageSource', 'logo-UBO-1-removebg-preview.png', 'Position', [150 650 100 50]);
            
            % Add Central Left LEMAR Logo (top-right corner)
            uiimage(app.UIFigure, 'ImageSource', 'logo-lemar-big.png', 'Position', [280 650 100 50]);

            % Main title
            uilabel(app.UIFigure, 'Text', 'Abalone DEB Simulation Framework v.1.3', ...
                'Position', [50 610 400 30], 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

            % Subtitle
            uilabel(app.UIFigure, 'Text', 'Authors: Paulo F. Lagos, Alva Curstodder, Laure Pecquerie', ...
                'Position', [50 598 400 20], 'FontSize', 12, 'HorizontalAlignment', 'center');
            uilabel(app.UIFigure, 'Text', 'Date: September 2024', 'Position', [50 585 400 20], 'FontSize', 10, 'HorizontalAlignment', 'center');

            % Create RunSimulationButton
            app.RunSimulationButton = uibutton(app.UIFigure, 'push', 'Text', 'Run Simulation', ...
                'Position', [100 10 200 40], 'ButtonPushedFcn', @(~,event)onRunSimulationButtonPressed(app, event), ...
                'BackgroundColor', [0.47 0.67 0.19], 'FontSize', 16, 'FontWeight', 'bold', 'FontColor', 'white');

            % Create InstructionsButton next to the RunSimulationButton
            app.InstructionsButton = uibutton(app.UIFigure, 'push', 'Text', 'Instructions', ...
                'Position', [320 10 150 40], 'ButtonPushedFcn', @(~, event)onInstructionsButtonPressed(app, event), ...
                'FontSize', 14, 'FontWeight', 'bold');
            
            % Create ScenarioListBox
            app.ScenarioListBox = uilistbox(app.UIFigure, 'Position', [100 489 150 70], ...
                'Items', {'Current', 'SSP1s19', 'SSP1s26', 'SSP1s45', 'SSP3sBaseline', 'SSP5sBaseline'}, ...
                'Multiselect', 'on', 'ValueChangedFcn', @(~,event)updateSelectedScenarios(app, event));  

            % Create LocationListBox
            app.LocationListBox = uilistbox(app.UIFigure, 'Position', [100 398 150 70], ...
                'Items', {'brest', 'roscoff', 'nantes', 'larochelle', 'med'}, ...
                'Multiselect', 'on', 'ValueChangedFcn', @(~,event)updateSelectedLocations(app, event));  

            % Create FoodQualityListBox
            app.FoodQualityListBox = uilistbox(app.UIFigure, 'Position', [100 307 150 70], ...
                'Items', {'0', '0.15','0.3'}, 'Multiselect', 'on', 'ValueChangedFcn', @(~,event)updateSelectedFoodQuality(app, event));  

            % Create EnvVariabilityListBox
            app.EnvVariabilityListBox = uilistbox(app.UIFigure, 'Position', [100 218 150 70], ...
                'Items', {'1', '2'}, 'Multiselect', 'on', 'ValueChangedFcn', @(~,event)updateSelectedEnvVariability(app, event));  
            
            % Create ZoomFactor Text Area for multiple inputs
            app.ZoomFactorListBox = uitextarea(app.UIFigure, 'Position', [100 175 150 22], 'ValueChangedFcn', @(~,event)updateSelectedZoomFactor(app, event));
            
            % Add note above the Instructions Button
            uilabel(app.UIFigure, 'Text', {'Note: To select multiple items', 'press Control (Ctrl) or Shift.'}, ...
                'Position', [320 60 150 40], 'FontSize', 10, 'FontColor', [0.5 0.5 0.5], 'BackgroundColor', app.UIFigure.Color);

            % Create TextAreas to display selected values
            app.SelectedScenariosLabel = uitextarea(app.UIFigure, 'Position', [260 489 200 70], 'Editable', 'off', 'Value', 'Selected Scenarios will appear here');
            app.SelectedLocationsLabel = uitextarea(app.UIFigure, 'Position', [260 398 200 70], 'Editable', 'off', 'Value', 'Selected Locations will appear here');
            app.SelectedFoodQualityLabel = uitextarea(app.UIFigure, 'Position', [260 307 200 70], 'Editable', 'off', 'Value', 'Selected Food Quality will appear here');
            app.SelectedEnvVariabilityLabel = uitextarea(app.UIFigure, 'Position', [260 218 200 70], 'Editable', 'off', 'Value', 'Selected Environmental Variability: (1) Constant or (2) Variable Environment');
            app.SelectedZoomFactorLabel  = uitextarea(app.UIFigure, 'Position', [260 175 200 22], 'Editable', 'off', 'Value', 'Input "," separate values');
            
            % Create labels and numeric fields
            uilabel(app.UIFigure, 'Text', 'Select climate scenario:', 'Position', [100 558 150 20]);
            uilabel(app.UIFigure, 'Text', 'Select locations:', 'Position', [100 465 150 20]);
            uilabel(app.UIFigure, 'Text', 'Enter food quality:', 'Position', [100 375 150 20]);
            uilabel(app.UIFigure, 'Text', 'Environmental variability:', 'Position', [100 285 150 20]);
            uilabel(app.UIFigure, 'Text', 'Enter number of individuals:', 'Position', [100 194 150 20]);
            uilabel(app.UIFigure, 'Text', 'Number of simulations in time:', 'Position', [100 155 200 20]);
            app.GenNumberNumericField = uieditfield(app.UIFigure, 'numeric', 'Position', [100 135 150 22], 'Value', 1);
            uilabel(app.UIFigure, 'Text', 'Simulation start day', 'Position', [100 110 200 20]);
            app.BornDateNumericField = uieditfield(app.UIFigure, 'numeric', 'Position', [100 90 150 22], 'Value', 152);
        end
    end

    % App initialization and construction
    methods (Access = public)

        % Construct app
        function app = AbaloneApp
            createComponents(app);
        end
    end
end
