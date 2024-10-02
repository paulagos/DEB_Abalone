classdef AbaloneApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        % The existing app UI components
        UIFigure                 matlab.ui.Figure
        RunSimulationButton      matlab.ui.control.Button
        ScenarioListBox          matlab.ui.control.ListBox
        LocationListBox          matlab.ui.control.ListBox
        FoodQualityListBox       matlab.ui.control.ListBox
        EnvVariabilityListBox    matlab.ui.control.ListBox
        SelectedScenariosLabel   matlab.ui.control.TextArea
        SelectedLocationsLabel   matlab.ui.control.TextArea
        SelectedFoodQualityLabel matlab.ui.control.TextArea
        SelectedEnvVariabilityLabel matlab.ui.control.TextArea
        ZoomFactorNumericField   matlab.ui.control.NumericEditField
        GenNumberNumericField    matlab.ui.control.NumericEditField
        BornDateNumericField     matlab.ui.control.NumericEditField
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

       
        % Step 1: Run Simulation Button callback
        function onRunSimulationButtonPressed(app, event)
            scenarios = app.ScenarioListBox.Value;
            locations = app.LocationListBox.Value;
            foodQuality = app.FoodQualityListBox.Value;
            envVariability = app.EnvVariabilityListBox.Value;
            zoomFactor = app.ZoomFactorNumericField.Value;
            genNumber = app.GenNumberNumericField.Value;
            bornDate = app.BornDateNumericField.Value;

            simu = app.prepareSimulationData(scenarios, locations, foodQuality, envVariability, zoomFactor, genNumber, bornDate);

            results = main(simu);

            app.displaySimulationResults(results);
        end

        % Step 2: Prepare Simulation Data
        function simu = prepareSimulationData(app, scenarios, locations, foodQuality, envVariability, zoomFactor, genNumber, bornDate)
            simu.scenarios = scenarios;
            simu.locations = locations;
            simu.foodQuality = foodQuality;
            simu.envVariability = envVariability;
            simu.zoomFactor = zoomFactor;
            simu.genNumber = genNumber;
            simu.bornDate = bornDate;

            simu.envT = 1;
            simu.envX = 1;
        end

        % Step 4: Display Simulation Results
        function displaySimulationResults(app, results)
            figure;
            plot(results.time, results.wetWeight);
            title('Simulation Results: Wet Weight Over Time');
            xlabel('Time (days)');
            ylabel('Wet Weight (g)');
            
         
        end
    end
   

    % Create UI components
    methods (Access = private)

        % Code to create app components
        function createComponents(app)

            % Create UIFigure and components
            app.UIFigure = uifigure('Position', [100 100 500 700]);

            % Main title
            uilabel(app.UIFigure, 'Text', 'Abalone DEB Simulation Framework v.1.0', ...
                'Position', [50 620 400 30], 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

            % Subtitle
            uilabel(app.UIFigure, 'Text', 'Authors: Dr. Paulo F. Lagos and Dr. Laure Pecquerie', ...
                'Position', [50 600 400 20], 'FontSize', 12, 'HorizontalAlignment', 'center');
            uilabel(app.UIFigure, 'Text', 'Date: September 2024', ...
                'Position', [50 585 400 20], 'FontSize', 12, 'HorizontalAlignment', 'center');

            % Create fun RunSimulationButton (Customized color)
            app.RunSimulationButton = uibutton(app.UIFigure, 'push', 'Text', 'Run Simulation', ...
                'Position', [100 10 200 40], 'ButtonPushedFcn', @(~,event)RunSimulationButtonPushed(app, event), ...
                'BackgroundColor', [0.47 0.67 0.19], 'FontSize', 14, 'FontWeight', 'bold', 'FontColor', 'white');

            % Adjust the position of labels for each option box (moved higher)
            uilabel(app.UIFigure, 'Text', 'Select Simulation Scenarios', ...
                'Position', [100 555 200 20]);

            uilabel(app.UIFigure, 'Text', 'Select Locations', ...
                'Position', [100 465 200 20]);

            uilabel(app.UIFigure, 'Text', 'Select Food Quality', ...
                'Position', [100 375 200 20]);

            uilabel(app.UIFigure, 'Text', 'Select Environmental Variability', ...
                'Position', [100 285 200 20]);

            % Create ScenarioListBox
            app.ScenarioListBox = uilistbox(app.UIFigure, 'Position', [100 500 150 70], ...
                'Items', {'Current', 'SSP1s19', 'SSP1s26', 'SSP1s45', 'SSP3sBaseline', 'SSP5sBaseline'}, ...
                'Multiselect', 'on', 'ValueChangedFcn', @(~,event)updateSelectedScenarios(app, event));  

            % Create LocationListBox
            app.LocationListBox = uilistbox(app.UIFigure, 'Position', [100 410 150 70], ...
                'Items', {'brest', 'roscoff', 'nantes', 'larochelle', 'med'}, ...
                'Multiselect', 'on', 'ValueChangedFcn', @(~,event)updateSelectedLocations(app, event));  

            % Create FoodQualityListBox
            app.FoodQualityListBox = uilistbox(app.UIFigure, 'Position', [100 320 150 70], ...
                'Items', {'0', '0.13'}, 'Multiselect', 'on', 'ValueChangedFcn', @(~,event)updateSelectedFoodQuality(app, event));  

            % Create EnvVariabilityListBox
            app.EnvVariabilityListBox = uilistbox(app.UIFigure, 'Position', [100 230 150 70], ...
                'Items', {'1', '2'}, 'Multiselect', 'on', 'ValueChangedFcn', @(~,event)updateSelectedEnvVariability(app, event));  

            % Create TextAreas to display selected values
            app.SelectedScenariosLabel = uitextarea(app.UIFigure, 'Position', [260 500 200 70], ...
                'Editable', 'off', 'Value', 'Selected Scenarios will appear here');
            app.SelectedLocationsLabel = uitextarea(app.UIFigure, 'Position', [260 410 200 70], ...
                'Editable', 'off', 'Value', 'Selected Locations will appear here');
            app.SelectedFoodQualityLabel = uitextarea(app.UIFigure, 'Position', [260 320 200 70], ...
                'Editable', 'off', 'Value', 'Selected Food Quality will appear here');
            app.SelectedEnvVariabilityLabel = uitextarea(app.UIFigure, 'Position', [260 230 200 70], ...
                'Editable', 'off', 'Value', 'Selected Environmental Variability will appear here');

            % Create labels and numeric fields for zoom, gen number, and born date
            uilabel(app.UIFigure, 'Text', 'Number of individuals', ...
                'Position', [100 205 200 20]);
            app.ZoomFactorNumericField = uieditfield(app.UIFigure, 'numeric', 'Position', [100 185 150 22], ...
                'Value', 1);

            uilabel(app.UIFigure, 'Text', 'Number of simulations in time', ...
                'Position', [100 165 200 20]);
            app.GenNumberNumericField = uieditfield(app.UIFigure, 'numeric', 'Position', [100 145 150 22], ...
                'Value', 1);

            uilabel(app.UIFigure, 'Text', 'Simulation start day', ...
                'Position', [100 125 200 20]);
            app.BornDateNumericField = uieditfield(app.UIFigure, 'numeric', 'Position', [100 105 150 22], ...
                'Value', 152);
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