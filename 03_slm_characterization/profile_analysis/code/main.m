% ----------------------------------------------------------
% Extract circular intensity profiles from lab and simulated images,
% fit profiles using Fourier model (fourier2), and export coefficients.
%
% Workflow:
%   1. Lab: Extract green channel, sample intensity on calibrated circle
%   2. Sim: User selects circle radius via ginput (saved to log file)
%   3. Fit: Fourier2 model fitting for all profiles
%   4. Compare: Generate overlay plots and metrics (AvgDiff, STD)
%   5. Export: Coefficients, figures, metrics, and LaTeX tables
%
% Input: ../preimages/ (lab images), ../Va_<va>_Vp_<vp>/ (sim images)
% Output: ../data/ (fit coefficients, metrics, LaTeX tables)
%         ../image/ (diagnostic figures, Fourier comparison plots)
%
% Notes: Lab circle radius/center are hard-coded (calibrated once).
%        Sim circle radius selected manually via ginput and saved to log file.
%        Generates LaTeX tables for AverageDifference and STD metrics.
%        Green highlighting indicates minimum absolute value per m.
% ----------------------------------------------------------

%% ------- [ Environment Setup ] -------
% Clear workspace, close figures, compact display
clc;
clear; close all; format compact; format bank;

%% ------- [ Working Directory Sync ] -------
% Auto-change current folder to the active script's location
disp('Auto change folder enabled.');
activeFile = matlab.desktop.editor.getActiveFilename;
if ~isempty(activeFile)
    cd(fileparts(activeFile));
    fprintf('Changed folder to: %s\n', pwd);
end

%% ------- [ Output Paths ] -------
% Define and create output directories if they do not exist
dirData = fullfile('..', 'data');
dirImage = fullfile('..', 'image');

% Create directories safely
if ~exist(dirData, 'dir')
    mkdir(dirData);
end
if ~exist(dirImage, 'dir')
    mkdir(dirImage);
end

%% ------- [ Parameters ] -------
% Case types for comparison
caseName = {'lab', 'sim'};

%% ------- [ Main Loop: Process m1 ] -------
% Loop over grating orders
m1Values = [10, 25, 50];

for m1 = 25%m1Values
    fprintf('Processing m1 = %d\n', m1);
    
    for caseNumber = 1:2
        switch caseName{caseNumber}
            
            case 'lab'
                %% ------- [ Lab Images ] -------
                % Read experimental image, extract green channel, normalize intensity
                
                name = sprintf('lab_m1_%d', m1);
                filePath = fullfile('..', 'preimages', [name '.png']);
                original_image = imread(filePath);        % Read RGB image
                image = double(original_image(:, :, 2));  % Use green channel
                image = image ./ max(image(:));           % Normalize to [0, 1]
                
                s = size(image);
                lx = s(2) / 2;
                ly = s(1) / 2;
                x0 = -lx:lx;
                y0 = -ly:ly;
                [x, y] = meshgrid(x0, y0);
                [tet, r] = cart2pol(x, y);
                
                %% ------- [ Circular Path: Lab ] -------
                % Calibrated radius and center (obtained once via manual selection)
                switch m1
                    case 10
                        radius = 89.6158;
                        MiddlePoint = [236.6353, 189.7653];
                    case 25
                        radius = 133.0573;
                        MiddlePoint = [236.2699, 191.7164];
                    case 50
                        radius = 183.7709;
                        MiddlePoint = [235.6723, 187.2927];
                end
                
                numPoints = max(s);  % Number of points on circular path
                
                % Generate circular path coordinates
                theta = linspace(0, 2*pi, numPoints);
                x_circle = radius * cos(theta) + MiddlePoint(1);
                y_circle = radius * sin(theta) + MiddlePoint(2);
                
                %% ------- [ Profile Extraction: Lab ] -------
                % Sample intensity along circular path
                fig = figure();
                s1 = subplot(2, 1, 1);
                s1.Units = "centimeters";
                
                imshow(original_image);
                hold on;
                plot(x_circle, y_circle, 'g', 'LineWidth', 1.5);
                set(gca, 'YDir', 'normal');
                
                profile1 = improfile(image, x_circle, y_circle);
                tet = linspace(-pi, pi, length(profile1));
                
                % Calculate visibility
                Imax = max(profile1);
                Imin = min(profile1);
                visibility = (Imax - Imin) / (Imax + Imin);
                fprintf('  Visibility: %.3f\n', visibility);
                
                s2 = subplot(2, 1, 2);
                s2.Units = "centimeters";
                s2.PlotBoxAspectRatio = s1.PlotBoxAspectRatio;
                xlabel('$x$ (px)', 'Interpreter', 'latex');
                ylabel('Intensity', 'Interpreter', 'latex');
                box on;
                
                %% ------- [ Fourier Fit: Lab ] -------
                % Fit profile using Fourier series model (fourier2)
                [fitResult, ~] = fit(tet', profile1, "fourier2");
                fprintf('  Fit coefficients saved.\n');
                
                % Save coefficients
                coeffData = coeffvalues(fitResult)';
                writematrix(coeffData, fullfile(dirData, [name '.txt']));
                
                % Plot fitted curve
                plot(tet, feval(fitResult, tet) ./ max(feval(fitResult, tet)), 'r-', 'LineWidth', 2);
                ylim([0, 1]);
                xlim([0, pi]);
                set(gca, 'defaultTextInterpreter', 'latex', 'TickLabelInterpreter', 'latex', 'FontSize', 15);
                xticks([0, pi]);
                xticklabels({'$0$', '$\pi$'});
                yticks([0, 1]);
                
                % Export figure
                exportgraphics(fig, fullfile(dirImage, [name '.pdf']));
                close(fig);
                
            case 'sim'
                %% ------- [ Sim Images ] -------
                % Loop over Va and Vp, read simulated intensity images
                
                Va = [0.1, 0.3, 0.5];
                Vp = [pi/4, pi/2, 3*pi/4];  % [0.7854, 1.5708, 2.3562]
                
                for va = Va
                    for vp = pi/4%Vp
                        fprintf('  Va = %.1f, Vp = %.2f pi\n', va, vp/pi);
                        
                        % Construct folder and file names using consistent convention
                        folderName = sprintf('Va_%.1f_Vp_%.2fpi', va, vp/pi);
                        folderPath = fullfile('../../generate_beam_profiles', folderName);
                        
                        % Read gray-scale image for analysis (saved as PNG from previous script)
                        fileName = sprintf('gray_propagation_m1=%d_Va=%.1f_Vp=%.2fpi.png', ...
                                           m1, va, vp/pi);
                        filePath = fullfile(folderPath, fileName);
                        
                        % If image is saved as PDF, convert to image first or save as PNG in previous script
                        % For now, assume PNG format
                        try
                            image = double(rgb2gray(imread(filePath)));
                            image = image ./ max(image(:));  % Normalize to [0, 1]
                        catch
                            warning('Could not read file: %s. Skipping...', filePath);
                            continue;
                        end
                        
                        % Also read hot-map image for reference (optional)
                        hotFileName = sprintf('hot_propagation_m1=%d_Va=%.1f_Vp=%.2fpi.png', ...
                                              m1, va, vp/pi);
                        hotFilePath = fullfile(folderPath, hotFileName);
                        if exist(hotFilePath, 'file')
                            original_image = imread(hotFilePath);
                        else
                            original_image = image;  % Use gray image as fallback
                        end
                        
                        s = size(image);
                        lx = s(2) / 2;
                        ly = s(1) / 2;
                        x0 = -lx:lx;
                        y0 = -ly:ly;
                        [x, y] = meshgrid(x0, y0);
                        [tet, r] = cart2pol(x, y);
                        
                        %% ------- [ Circular Path: Sim ] -------
                        % Radius selection with option to use existing or select new
                        
                        % Check if log file exists
                        radiusLogFile = fullfile(dirData, 'radius_calibration_log.txt');
                        logExists = exist(radiusLogFile, 'file');
                        
                        % Ask user preference
                        if logExists
                            % Check if radius exists for current combination in log file
                            % Read the log file and search for current combination
                            fid = fopen(radiusLogFile, 'r');
                            logContent = textscan(fid, '%s', 'Delimiter', '\n');
                            fclose(fid);
                            logLines = logContent{1};
                            
                            % Search for current combination
                            searchPattern = sprintf('m1=%d, Va=%.1f, Vp=%.2fpi', m1, va, vp/pi);
                            foundRadius = false;
                            existingRadius = 0;
                            
                            for i = 1:length(logLines)
                                if contains(logLines{i}, searchPattern)
                                    % Extract radius value from the line
                                    parts = strsplit(logLines{i}, '-> radius = ');
                                    if length(parts) == 2
                                        existingRadius = str2double(parts{2});
                                        foundRadius = true;
                                        break;
                                    end
                                end
                            end
                            
                            if foundRadius
                                choice = questdlg(sprintf('Radius found for m1=%d, Va=%.1f, Vp=%.2fpi:\nExisting radius = %.4f\n\nDo you want to use this radius or select a new one?', ...
                                                   m1, va, vp/pi, existingRadius), ...
                                                  'Radius Selection', ...
                                                  'Use Existing', 'Select New', 'Cancel', ...
                                                  'Use Existing');
                            else
                                choice = questdlg(sprintf('No radius found for m1=%d, Va=%.1f, Vp=%.2fpi\n\nDo you want to select a new radius?', ...
                                                   m1, va, vp/pi), ...
                                                  'Radius Selection', ...
                                                  'Select New', 'Cancel', ...
                                                  'Select New');
                            end
                        else
                            % No log file exists, must select new
                            choice = 'Select New';
                            fprintf('  No radius log file found. Selecting new radius...\n');
                        end
                        
                        % Process user choice
                        switch choice
                            case 'Use Existing'
                                % Use existing radius from log file
                                radius = existingRadius;
                                fprintf('  Using existing radius: %.4f\n', radius);
                                
                            case 'Select New'
                                % Display image and wait for user click on the ring
                                figRadius = figure('Name', sprintf('m1=%d, Va=%.1f, Vp=%.2fpi - Click on ring', m1, va, vp/pi), ...
                                                   'NumberTitle', 'off');
                                imagesc(image);
                                colormap hot;
                                axis image;
                                title(sprintf('m1=%d, Va=%.1f, Vp=%.2fpi\nClick on the ring', m1, va, vp/pi), ...
                                      'FontSize', 12);
                                hold on;
                                
                                % Mark image center
                                point1 = [s(2)/2, s(1)/2];  % Image center
                                plot(point1(1), point1(2), 'b+', 'MarkerSize', 10, 'LineWidth', 2);
                                
                                % User clicks on the ring
                                [clickX, clickY] = ginput(1);
                                point2 = [clickX, clickY];
                                
                                % Calculate radius
                                radius = sqrt((point2(1) - point1(1))^2 + (point2(2) - point1(2))^2);
                                
                                % Draw the circle on the image for verification
                                theta_circle = linspace(0, 2*pi, 100);
                                x_circle_verify = radius * cos(theta_circle) + point1(1);
                                y_circle_verify = radius * sin(theta_circle) + point1(2);
                                plot(x_circle_verify, y_circle_verify, 'g-', 'LineWidth', 2);
                                plot(point2(1), point2(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
                                
                                % Print the radius value for use in getRadius function
                                fprintf('  >>> NEW RADIUS: m1=%d, Va=%.1f, Vp=%.2fpi -> radius = %.4f\n', ...
                                        m1, va, vp/pi, radius);
                                
                                % Check if this combination already exists in log file and update or append
                                if logExists
                                    % Read existing log file
                                    fid = fopen(radiusLogFile, 'r');
                                    logContent = textscan(fid, '%s', 'Delimiter', '\n');
                                    fclose(fid);
                                    logLines = logContent{1};
                                    
                                    % Check if current combination exists
                                    searchPattern = sprintf('m1=%d, Va=%.1f, Vp=%.2fpi', m1, va, vp/pi);
                                    foundIndex = -1;
                                    for i = 1:length(logLines)
                                        if contains(logLines{i}, searchPattern)
                                            foundIndex = i;
                                            break;
                                        end
                                    end
                                    
                                    if foundIndex > 0
                                        % Update existing entry
                                        logLines{foundIndex} = sprintf('m1=%d, Va=%.1f, Vp=%.2fpi -> radius = %.4f', ...
                                                                       m1, va, vp/pi, radius);
                                        % Write updated log file
                                        fid = fopen(radiusLogFile, 'w');
                                        for i = 1:length(logLines)
                                            fprintf(fid, '%s\n', logLines{i});
                                        end
                                        fclose(fid);
                                        fprintf('  Updated existing entry in log file.\n');
                                    else
                                        % Append new entry
                                        fid = fopen(radiusLogFile, 'a');
                                        fprintf(fid, 'm1=%d, Va=%.1f, Vp=%.2fpi -> radius = %.4f\n', ...
                                                m1, va, vp/pi, radius);
                                        fclose(fid);
                                        fprintf('  Added new entry to log file.\n');
                                    end
                                else
                                    % Create new log file with first entry
                                    fid = fopen(radiusLogFile, 'w');
                                    fprintf(fid, 'm1=%d, Va=%.1f, Vp=%.2fpi -> radius = %.4f\n', ...
                                            m1, va, vp/pi, radius);
                                    fclose(fid);
                                    fprintf('  Created new log file with first entry.\n');
                                end
                                
                                % Wait for user to confirm before continuing
                                disp('  Press any key to continue...');
                                pause;
                                
                                % Close the figure
                                close(figRadius);
                                
                            case 'Cancel'
                                % User cancelled, skip this iteration
                                fprintf('  Skipping m1=%d, Va=%.1f, Vp=%.2fpi\n', m1, va, vp/pi);
                                continue;
                                
                            otherwise
                                % Default: skip
                                fprintf('  Skipping m1=%d, Va=%.1f, Vp=%.2fpi\n', m1, va, vp/pi);
                                continue;
                        end
                        
                        numPoints = max(s);  % Number of points on circular path
                        
                        % Generate circular path coordinates centered at image center
                        theta = linspace(0, 2*pi, numPoints);
                        x_circle = radius * cos(theta) + lx;
                        y_circle = radius * sin(theta) + ly;
                        
                        %% ------- [ Profile Extraction: Sim ] -------
                        % Create figure with image and profile
                        fig = figure();
                        s1 = subplot(2, 1, 1);
                        s1.Units = "centimeters";
                        
                        imagesc(image);
                        colormap hot;
                        xticks([]); yticks([]);
                        axis image;
                        hold on;
                        plot(x_circle, y_circle, 'g', 'LineWidth', 1.5);
                        set(gca, 'YDir', 'normal');
                        
                        % Sample intensity along circular path
                        profile1 = improfile(image, x_circle, y_circle);
                        tet = linspace(-pi, pi, length(profile1));
                        
                        s2 = subplot(2, 1, 2);
                        s2.Units = "centimeters";
                        
                        %% ------- [ Fourier Fit: Sim ] -------
                        % Fit using Fourier series model
                        [fitResult, ~] = fit(tet', profile1, "fourier2");
                        
                        % Save coefficients
                        simName = sprintf('sim_m1_%d_Va_%.1f_Vp_%.2fpi', m1, va, vp/pi);
                        coeffData = coeffvalues(fitResult)';
                        writematrix(coeffData, fullfile(dirData, [simName '.txt']));
                        
                        % Plot fitted curve
                        plot(tet, feval(fitResult, tet) ./ max(feval(fitResult, tet)), 'r-', 'LineWidth', 2);
                        xlim([0, pi]);
                        ylim([0, 1]);
                        set(gca, 'defaultTextInterpreter', 'latex', 'TickLabelInterpreter', 'latex', 'FontSize', 15);
                        xticks([0, pi]);
                        xticklabels({'$0$', '$\pi$'});
                        yticks([0, 1]);
                        xlabel('$x$ (px)', 'Interpreter', 'latex');
                        ylabel('Intensity', 'Interpreter', 'latex');
                        box on;
                        s2.PlotBoxAspectRatio = s1.PlotBoxAspectRatio;
                        
                        % Export figure
                        exportgraphics(fig, fullfile(dirImage, [simName '.pdf']));
                        close(fig);
                    end
                end
        end
    end
end


%% ------- [ Comparison: Fourier2 Profiles ] -------
% Compare reconstructed Fourier2 profiles for lab vs sim cases

% Parameters: Angular sampling
% Number of angular samples for profile reconstruction
N = 2^10;                            % Number of angular samples
tet = linspace(-pi, pi, N);          % Angular coordinate (rad)

% Visualization settings
% Case labels and plot colors
caseName = {'lab', 'sim'};

customColors = [
    1.00, 0.00, 0.00;  % Red
    0.00, 1.00, 0.00;  % Green
    0.00, 0.00, 1.00;  % Blue
    1.00, 0.65, 0.00;  % Orange
    0.75, 0.75, 0.75;  % Light Gray
    0.60, 0.00, 0.60;  % Purple
    0.00, 0.75, 1.00;  % Cyan
    0.60, 0.60, 0.00;  % Olive
    0.50, 0.00, 0.50;  % Dark Purple
    0.75, 0.00, 0.75;  % Pink
];

% Main comparison loop
% Compare reconstructed Fourier2 profiles for each m1
for m1 = m1Values
    fprintf('Comparing Fourier profiles for m1 = %d\n', m1);
    
    fig = figure('Position', [418, 132, 992, 520]);
    hold on;
    colorIndex = 1;  % Start with the first color in the array
    
    % Initialize result cell array for all curves
    result = cell(10, 1);  % 1 lab + 9 sim curves
    
    for caseNumber = 1:2
        switch caseName{caseNumber}
            case 'lab'
                % Lab profile: Read Fourier2 coefficients, reconstruct, normalize, and plot
                filename = sprintf('lab_m1_%d.txt', m1);
                filePath = fullfile(dirData, filename);
                
                if exist(filePath, 'file')
                    data = ExtractCoefficient(filePath);
                    result{1} = Fourier2nd(tet, data);
                    result{1} = result{1} ./ max(result{1});
                    plot(tet, result{1}, 'k', 'LineWidth', 2, 'DisplayName', 'Lab. data');
                    fprintf('  Loaded lab data for m1 = %d\n', m1);
                else
                    warning('Lab data not found for m1 = %d: %s', m1, filePath);
                end
                
            case 'sim'
                % Sim profiles: Loop over (Va, Vp), reconstruct 9 profiles, normalize, and plot
                counter = 1;
                
                Va_list = [0.1, 0.3, 0.5];
                Vp_list = [pi/4, pi/2, 3*pi/4];  % [0.7854, 1.5708, 2.3562]
                
                for va = Va_list
                    for vp = Vp_list
                        simName = sprintf('sim_m1_%d_Va_%.1f_Vp_%.2fpi', m1, va, vp/pi);
                        filePath = fullfile(dirData, [simName '.txt']);
                        
                        if exist(filePath, 'file')
                            data = ExtractCoefficient(filePath);
                            result{counter + 1} = Fourier2nd(tet, data);
                            result{counter + 1} = result{counter + 1} ./ max(result{counter + 1});
                            
                            plot(tet, result{counter + 1}, ...
                                'Color', customColors(colorIndex, :), ...
                                'LineStyle', '-', ...
                                'LineWidth', 1, ...
                                'DisplayName', sprintf('$V_a = %.1f$, $V_p = %.2f\\pi$', va, vp/pi));
                            
                            colorIndex = colorIndex + 1;  % Move to the next color
                            counter = counter + 1;
                        else
                            warning('Sim data not found: %s', filePath);
                        end
                    end
                end
        end
    end
    
    hold off;
    
    % Plot formatting: Axes, legend, ticks, and export
    ylabel('Normalized intensity', 'Interpreter', 'latex');
    xlabel('$\theta$ (rad)', 'Interpreter', 'latex');
    legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'best');
    set(gca, 'defaulttextinterpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gca, 'fontsize', 18);
    ylim([0, 1.1]);
    
    % Set x-axis limits based on the lab data's angular frequency
    if exist('data', 'var') && ~isempty(data)
        w = data(6);
        xLim = 2*pi / w;
        xlim([-xLim - 0.05, xLim + 0.05]);
        
        % Set x-axis ticks in terms of pi
        tickVal = round(xLim, 2);
        xticks(sort([-tickVal, 0, tickVal]));
        xticklabels({[num2str(tickVal) '$\pi$'], '$0$', [num2str(-tickVal) '$\pi$']});
    else
        xlim([-pi, pi]);
        xticks([-pi, 0, pi]);
        xticklabels({'$-\pi$', '$0$', '$\pi$'});
    end
    
    box on;
    
    % Export figure
    fileName = sprintf('Fourier_comparison_m1_%d.pdf', m1);
    exportgraphics(fig, fullfile(dirImage, fileName));
    fprintf('  Exported: %s\n', fileName);
    close(fig);
    
    % Metrics: Difference statistics
    % Compute difference between sim profiles and lab profile
    
    % Check if lab data exists
    if exist('result', 'var') && ~isempty(result{1})
        % Legacy metrics (unlabeled)
        AverageDifference = zeros(9, 1);
        STD = zeros(9, 1);
        
        for jj = 1:9
            if ~isempty(result{jj + 1})
                diff = result{jj + 1} - result{1};
                AverageDifference(jj) = mean(diff);
                STD(jj) = std(diff);
            end
        end
        
        AverageDifference = round(AverageDifference, 4);
        STD = round(STD, 4);
        
        % Save legacy vectors
        writematrix(AverageDifference, fullfile(dirData, sprintf('m1_%d_AverageDifference.txt', m1)));
        writematrix(STD, fullfile(dirData, sprintf('m1_%d_STD.txt', m1)));
        
        % Metrics: Labeled table
        % Difference statistics versus lab (save as labeled table)
        Va_list = [0.1, 0.3, 0.5];
        Vp_list = [pi/4, pi/2, 3*pi/4];
        
        row = 0;
        Va_col   = zeros(numel(Va_list) * numel(Vp_list), 1);
        Vp_col   = zeros(numel(Va_list) * numel(Vp_list), 1);
        AvgDiff  = zeros(numel(Va_list) * numel(Vp_list), 1);
        StdDiff  = zeros(numel(Va_list) * numel(Vp_list), 1);
        
        for ia = 1:numel(Va_list)
            for ip = 1:numel(Vp_list)
                row = row + 1;
                simIdx = (ia - 1) * numel(Vp_list) + ip + 1;  % result{2}..result{10}
                
                if ~isempty(result{simIdx})
                    diff = result{simIdx} - result{1};
                    Va_col(row)  = Va_list(ia);
                    Vp_col(row)  = Vp_list(ip);
                    AvgDiff(row) = mean(diff);
                    StdDiff(row) = std(diff);
                end
            end
        end
        
        AvgDiff = round(AvgDiff, 4);
        StdDiff = round(StdDiff, 4);
        
        T = table(Va_col, Vp_col, AvgDiff, StdDiff, ...
            'VariableNames', {'Va', 'Vp', 'AverageDifference', 'STD'});
        % Save as tab-delimited text file
        writetable(T, fullfile(dirData, sprintf('m1_%d_SummaryMetrics.txt', m1)), 'Delimiter', '\t');
        
        % Save as Excel file
        writetable(T, fullfile(dirData, sprintf('m1_%d_SummaryMetrics.xlsx', m1)));

        fprintf('  Saved summary metrics for m1 = %d\n', m1);
    else
        warning('Lab data not available for m1 = %d. Skipping metrics.', m1);
    end
end

%% ------- [ Export LaTeX Tables from Individual Files ] -------

% Vp labels in the exact format you want
Vp_labels = {'$\nicefrac{\pi}{4}$', '$\nicefrac{\pi}{2}$', '$\nicefrac{3\pi}{4}$'};
Vp_labels_std = {'$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$'};

fid = fopen(fullfile(dirData, 'SummaryMetrics_All.tex'), 'w');

% ============================================================
% TABLE 1: Average Difference
% ============================================================
fprintf(fid, '\\begin{table}[!b]\n');
fprintf(fid, '\t\\caption{\\label{tab1}Mean differences between simulated and experimental intensity patterns for various values of $m$, $V_a$, and $V_p$.}\n');
fprintf(fid, '\t\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lcccc}\n');
fprintf(fid, '\t\t\\toprule\n');
fprintf(fid, '\t\t$m$ & \\diagbox{$V_p$}{$V_a$} & $0.1$ & $0.3$ & $0.5$ \\\\\n');
fprintf(fid, '\t\t\\midrule\n');

% Loop over m1 values
for mIdx = 1:length(m1Values)
    m1 = m1Values(mIdx);
    
    % Read the AverageDifference file
    avgFile = fullfile(dirData, sprintf('m1_%d_AverageDifference.txt', m1));
    if exist(avgFile, 'file')
        AvgDiff = readmatrix(avgFile);
        % Find minimum ABSOLUTE value for this m value
        [~, minAbsIdx] = min(abs(AvgDiff));
        minAvgDiff = AvgDiff(minAbsIdx);  % Actual value at the minimum absolute position
    else
        warning('File not found: %s', avgFile);
        continue;
    end
    
    fprintf(fid, '\t\t\\multirow{3}{*}{$%d$} \n', m1);
    
    % Order: Va outer loop, Vp inner loop
    for ip = 1:3  % Vp loop (rows)
        fprintf(fid, '\t\t& %s ', Vp_labels{ip});
        
        for ia = 1:3  % Va loop (columns)
            idx = (ia - 1) * 3 + ip;
            val = AvgDiff(idx);
            
            % Highlight if this is the minimum ABSOLUTE value for this m value
            if val == minAvgDiff
                fprintf(fid, '& \\colorbox{green!25}{$%.4f$} ', val);
            else
                fprintf(fid, '& $%.4f$ ', val);
            end
        end
        fprintf(fid, '\\\\\n');
    end
    
    if mIdx < length(m1Values)
        fprintf(fid, '\t\t\\midrule\n');
    end
end

fprintf(fid, '\t\t\\bottomrule\n');
fprintf(fid, '\t\\end{tabular*}\n');
fprintf(fid, '\\end{table}\n\n');

% ============================================================
% TABLE 2: STD
% ============================================================
fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\t\\caption{\\label{tab2}STDs of differences between simulated and experimental intensity patterns for various values of $m$, $V_a$, and $V_p$.}\n');
fprintf(fid, '\t\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lcccc}\n');
fprintf(fid, '\t\t\\toprule\n');
fprintf(fid, '\t\t$m$ & \\diagbox{$V_p$}{$V_a$} & $0.1$ & $0.3$ & $0.5$ \\\\\n');
fprintf(fid, '\t\t\\midrule\n');

for mIdx = 1:length(m1Values)
    m1 = m1Values(mIdx);
    
    stdFile = fullfile(dirData, sprintf('m1_%d_STD.txt', m1));
    if exist(stdFile, 'file')
        StdDiff = readmatrix(stdFile);
        % Find minimum ABSOLUTE value for this m value
        [~, minAbsIdx] = min(abs(StdDiff));
        minStd = StdDiff(minAbsIdx);  % Actual value at the minimum absolute position
    else
        warning('File not found: %s', stdFile);
        continue;
    end
    
    fprintf(fid, '\t\t\\multirow{3}{*}{$%d$} \n', m1);
    
    for ip = 1:3
        fprintf(fid, '\t\t& %s ', Vp_labels_std{ip});
        
        for ia = 1:3
            idx = (ia - 1) * 3 + ip;
            val = StdDiff(idx);
            
            % Highlight if this is the minimum ABSOLUTE value for this m value
            if val == minStd
                fprintf(fid, '& \\colorbox{green!25}{$%.4f$} ', val);
            else
                fprintf(fid, '& $%.4f$ ', val);
            end
        end
        fprintf(fid, '\\\\\n');
    end
    
    if mIdx < length(m1Values)
        fprintf(fid, '\t\t\\midrule\n');
    end
end

fprintf(fid, '\t\t\\bottomrule\n');
fprintf(fid, '\t\\end{tabular*}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid);
fprintf('Saved LaTeX tables to: %s\n', fullfile(dirData, 'SummaryMetrics_All.tex'));

%% ------- [ Local Functions ] -------

function data = ExtractCoefficient(filename)
    % Read Fourier2 coefficients stored as a 6-element vector: [a0 a1 b1 a2 b2 w]
    data = importdata(filename);
end

function result = Fourier2nd(tet, data)
    % Reconstruct the Fourier2 model:
    % f(tet) = a0 + a1*cos(w*tet) + b1*sin(w*tet) + a2*cos(2*w*tet) + b2*sin(2*w*tet)
    
    a0 = data(1);
    a1 = data(2);
    b1 = data(3);
    a2 = data(4);
    b2 = data(5);
    w = data(6);
    
    result = a0 + a1 * cos(tet * w) + b1 * sin(tet * w) + ...
             a2 * cos(2 * tet * w) + b2 * sin(2 * tet * w);
end