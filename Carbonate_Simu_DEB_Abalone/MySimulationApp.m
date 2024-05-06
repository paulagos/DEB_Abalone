classdef MySimulationApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure          matlab.ui.Figure
        RunSimulationButton  matlab.ui.control.Button
        ScenarioDropDown  matlab.ui.control.DropDown
        LocationDropDown  matlab.ui.control.DropDown
        FoodQualityDropDown  matlab.ui.control.DropDown
        EnvVariabilityDropDown  matlab.ui.control.DropDown
        ZoomFactorNumericField matlab.ui.control.NumericEditField
        StartYearNumericField matlab.ui.control.NumericEditField
        EndYearNumericField matlab.ui.control.NumericEditField
        
        % Labels as properties
        ScenarioDropDownLabel matlab.ui.control.Label
        LocationDropDownLabel matlab.ui.control.Label
        FoodQualityDropDownLabel matlab.ui.control.Label
        EnvVariabilityDropDownLabel matlab.ui.control.Label
        ZoomFactorLabel matlab.ui.control.Label 
        StartYearLabel matlab.ui.control.Label
        EndYearLabel matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: RunSimulationButton
        function RunSimulationButtonPushed(app, event)
            scenario = app.ScenarioDropDown.Value;
            location = app.LocationDropDown.Value;
            foodQuality = str2double(app.FoodQualityDropDown.Value);
            envVar = app.EnvVariabilityDropDown.Value;
            zoomFactor = app.ZoomFactorNumericField.Value;
            startYear = app.StartYearNumericField.Value;
            endYear = app.EndYearNumericField.Value;
            
            % Validate the input years
            if endYear <= startYear
                uialert(app.UIFigure, 'End Year must be greater than Start Year.', 'Input Error');
                return;
            end
            
            % Now call the main function with these parameters
            main(scenario, location, foodQuality, envVar, zoomFactor, genNumber);
            
            disp('Simulation completed');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)
            
            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'Abalone Simulation App';
            
            % LABELS for dropDown Menus
            % Create Scenario Dropdown Label
            app.ScenarioDropDownLabel = uilabel(app.UIFigure);
            app.ScenarioDropDownLabel.Position = [100 420 120 25];
            app.ScenarioDropDownLabel.Text = 'SSPs Climate Scenario:';
            
            % Create Locations Dropdown Label
            app.LocationDropDownLabel = uilabel(app.UIFigure);
            app.LocationDropDownLabel.Position = [100 370 120 22];
            app.LocationDropDownLabel.Text = 'Locations:';
                        
            % Create FoodQuality Dropdown Label
            app.FoodQualityDropDownLabel = uilabel(app.UIFigure);
            app.FoodQualityDropDownLabel.Position = [100 320 120 22];
            app.FoodQualityDropDownLabel.Text = 'Drop in Food quality (%):';
            
            % Create Environment Dropdow Label
            app.EnvVariabilityDropDownLabel = uilabel(app.UIFigure);
            app.EnvVariabilityDropDownLabel.Position = [100 270 120 22];
            app.EnvVariabilityDropDownLabel.Text = 'Environment type:';
            
            % StartYear Lable
            app.StartYearLabel = uilabel(app.UIFigure);
            app.StartYearLabel.Position = [100, 180, 100, 22];
            app.StartYearLabel.Text = 'Start Year:';
            
            % EndYear Label
            app.EndYearLabel = uilabel(app.UIFigure);
            app.EndYearLabel.Position = [220, 180, 100, 22];
            app.EndYearLabel.Text = 'End Year:';
            
            %%% DropDown Menus %%%
            % Create ScenarioDropDown 
            app.ScenarioDropDown = uidropdown(app.UIFigure);
            app.ScenarioDropDown.Position = [100 400 120 22];
            app.ScenarioDropDown.Items = {'Current', 'SSP1 1.9', 'SSP1 2.6','SSP1 4.5', 'SSP3 Baseline', 'SSP5 Baseline', 'All Scenarios'};
            
            % Create LocationsDropDown 
            app.LocationDropDown = uidropdown(app.UIFigure);
            app.LocationDropDown.Position = [100 350 120 22];
            app.LocationDropDown.Items = {'Bay of Brest', 'Bay of Morlaix', 'Bay of Villiane','La Rochelle','All locations'};
            
            % Create FoodQualityDropDown 
            app.FoodQualityDropDown = uidropdown(app.UIFigure);
            app.FoodQualityDropDown.Position = [100, 300, 120, 22];
            app.FoodQualityDropDown.Items = {'0', '15', '30'}; 
            
            % Create EnvVariabilityDropDown
            app.EnvVariabilityDropDown = uidropdown(app.UIFigure);
            app.EnvVariabilityDropDown.Position = [100, 250, 120, 22];
            app.EnvVariabilityDropDown.Items = {'Constant', 'Variable'};
            
            % Create StartYearNumericField  
            app.StartYearNumericField = uieditfield(app.UIFigure, 'numeric');
            app.StartYearNumericField.Position = [160 180 50 22];
            app.StartYearNumericField.Limits = [0 Inf];
            app.StartYearNumericField.Value = 1;

            % Create EndYearNumericField 
            app.EndYearNumericField = uieditfield(app.UIFigure, 'numeric');
            app.EndYearNumericField.Position = [280 180 50 22];
            app.EndYearNumericField.Limits = [0 Inf];
            app.EndYearNumericField.Value = 2;
            
            % Repeat this for other dropdowns, positioning them accordingly

            % Create RunSimulationButton
            app.RunSimulationButton = uibutton(app.UIFigure, 'push');
            app.RunSimulationButton.Position = [100 50 100 22];
            app.RunSimulationButton.Text = 'Run Simulation';
            app.RunSimulationButton.ButtonPushedFcn = createCallbackFcn(app, @RunSimulationButtonPushed, true);

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MySimulationApp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end