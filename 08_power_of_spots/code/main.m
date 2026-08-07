% ----------------------------------------------------------
% Simulate two-radial-grating diffraction and propagation through lens.
% Quantify central-spot metrics for m1 = [1,3,5]*m2 over m2 sweep.
%
% Workflow: Build gratings -> Propagate to second grating ->
%           Propagate through lens -> Detect central spot ->
%           Extract metrics (total power, density, intensity, area) ->
%           Compare simulation with experimental lab images
%
% Output: ../image/ (per-m2 bar plots, grouped bar charts, line plots, sim vs lab comparison)
%         ../data/  (central spot matrices)
%         ../result/ (text reports, aggregate tables, lab metrics)
%
% Notes: m2 = 10, m1 = [1,3,5]*m2 = [10, 30, 50].
%        Normalization relative to m1 = m2 (m1 = 10).
%        Lab images in ../labImage/ with naming m1=#.png.
% ----------------------------------------------------------

%% ------- [ Environment Setup ] -------
% Clear workspace, close figures, compact display

clc;
clear; close all; format compact; format bank;

%% ------- [ Output Paths ] -------
% Define and create output directories

dirImage  = fullfile('..', 'image');
dirData   = fullfile('..', 'data');
dirResult = fullfile('..', 'result');

% Create directories safely
if ~exist(dirImage, 'dir'),  mkdir(dirImage);  end
if ~exist(dirData, 'dir'),   mkdir(dirData);   end
if ~exist(dirResult, 'dir'), mkdir(dirResult); end

%% ------- [ Spatial Coordinates ] -------
% Build Cartesian grid and corresponding polar coordinates

Nx = 6000 + 1;      % Grid size along x
Ny = 6000 + 1;      % Grid size along y

lx = 1.5e-2;        % Half-window size along x [m]
ly = 1.2e-2;        % Half-window size along y [m]

x0 = linspace(-lx, lx, Nx);
y0 = linspace(-ly, ly, Ny);
[x, y] = meshgrid(x0, y0);
[tet, r] = cart2pol(x, y);  % tet: azimuth [rad], r: radius [m]

%% ------- [ Frequency Coordinates ] -------
% Frequency axes for angular-spectrum propagation

Lx = x0(end) - x0(1);   % Window width along x [m]
Ly = y0(end) - y0(1);   % Window height along y [m]
dx = Lx / Nx;           % Sample spacing along x [m]
dy = Ly / Ny;           % Sample spacing along y [m]

dfx = 1 / (Nx * dx);    % Fundamental frequency along x [1/m]
dfy = 1 / (Ny * dy);    % Fundamental frequency along y [1/m]

fx0 = dfx .* (-Nx/2 : Nx/2 - 1);
fy0 = dfy .* (-Ny/2 : Ny/2 - 1);
fx0 = fftshift(fx0);
fy0 = fftshift(fy0);
[fx, fy] = meshgrid(fx0, fy0);

%% ------- [ Input Field: SLM Aperture ] -------
% Rectangular SLM clear aperture with plane wave illumination

uPlane = zeros(size(x));
wx = 1.6e-2;        % SLM width [m]
wy = 1.2e-2;        % SLM height [m]
uPlane(x > -wx/2 & x < wx/2 & y > -wy/2 & y < wy/2) = 1;

%% ------- [ Optical Constants ] -------
% Define complex unit, wavelength, and wavenumber

i = complex(0, 1);
lambda = 532e-9;    % Wavelength in meters (532 nm)
k0 = 2*pi/lambda;   % Wavenumber [1/m]

%% ------- [ Grating Parameters ] -------
% Modulation parameters for hybrid grating

Vp = pi/2;          % Phase modulation amplitude [rad]
Va = 0.1;           % Amplitude visibility

%% ------- [ Main Simulation Loop ] -------
% Iterate over m2 values and analyze central spot metrics

m2Values = 10;
allResults = [];  % [m2, m1, RelTotalPower, RelPowerDensity, RelOriginIntensity, RelArea]

for m2Idx = 1:length(m2Values)
    
    m2 = m2Values(m2Idx);
    fprintf('m2 = %d\n', m2);

    m1Values = m2 * [1 3 5];
    analysis = struct();

    for m1Idx = 1:length(m1Values)
    
        m1 = m1Values(m1Idx);
        
        % Create first hybrid radial grating (amplitude + phase)
        t1 = (0.5 * (1 + Va * sign(cos(m1 * tet)))) .* ...
             exp(1i * Vp * sign(cos(m1 * tet)));
        
        % Propagate from SLM to second grating plane (z = 2.3 m)
        zToSecondGrating = 2.3;  % [m]
        uAtSecondGrating = ifft2(fft2(t1 .* uPlane) .* exp(-1i * pi * lambda * zToSecondGrating * (fx.^2 + fy.^2)));
     
        % Apply second binary radial amplitude grating (m2)
        t2 = 0.5 * (1 + sign(cos(m2 * (tet))));
        uAfterSecondGrating = uAtSecondGrating .* t2;
        
        % Propagate from second grating to lens plane (z = 5 cm)
        zToLens = 5e-2;  % [m]
        uAtLens = ifft2(fft2(uAfterSecondGrating) .* exp(-1i * pi * lambda * zToLens * (fx.^2 + fy.^2)));
        
        % Apply thin lens phase
        focalLength = 20e-2;  % Focal length [m]
        lensPhase = exp(-1i * k0 * r.^2 / (2 * focalLength));
        uAfterLens = uAtLens .* lensPhase;
        
        % Propagate from lens to detection plane (z = f + 1 cm)
        zDetection = 21e-2;  % [m]
        uDetection = ifft2(fft2(uAfterLens) .* exp(-1i * pi * lambda * zDetection * (fx.^2 + fy.^2)));
        
        % Compute intensity at detection plane
        intensity = abs(uDetection).^2;
        
        % ------- [ Display Intensity Pattern ] -------
        fig = figure('Visible', 'on', 'NumberTitle', 'off', 'Name', ['Figure ' num2str(m1Idx)]);
        set(fig, 'Position', [85.00 29.00 833.00 782.00]);

        s1 = subplot(4,3,1);
        imagesc(x0*1e+3, y0*1e+3, intensity);
        axis image;
        set(gca, 'YDir', 'normal'); colormap(hot); colorbar();
        sgtitle(['$m_1 = ', num2str(m1), '$, $m_2 = ', num2str(m2), ...
               '$, $z_{\mathrm{Detection}} = ', num2str(zDetection * 100), ...
               '$ cm, $f = ', num2str(focalLength * 100), '$ cm'], 'Interpreter', 'latex');
        xLim = 0.4; % [mm]
        yLim = 0.3; % [mm]
        xlim([-xLim, xLim]); xticks(-xLim:xLim:xLim);
        ylim([-yLim, yLim]); yticks(-yLim:yLim:yLim);
        clim([0 Inf])
        xlabel('$x$ (mm)');
        ylabel('$y$ (mm)');
        
        % ------- [ Analyze Central Spot Intensity and Area ] -------
        % Draws dashed guide lines on the central vertical and horizontal axes

        hold on;
        
        % Plot vertical dashed line at x = 0
        xCenter = x0(round(Nx/2));
        plot([xCenter xCenter]*1e3, [-yLim yLim], 'w--', 'LineWidth', 1.2);
        
        % Plot horizontal dashed line at y = 0
        yCenter = y0(round(Ny/2));
        plot([-xLim xLim], [yCenter yCenter]*1e3, 'w--', 'LineWidth', 1.2);
        
        % Extract intensity profile along vertical centerline (x = 0)
        intensityCol = intensity(:, round(Nx/2));
        
        % Plot vertical profile in subplot (top-right)
        s2 = subplot(4,3,2);
        plot(intensityCol, y0*1e3, 'k', 'LineWidth', 1.2);
        ylim([-yLim yLim]);
        xlabel('Intensity along $x = 0$');
        ylabel('$y$ (mm)');
        yticks(-yLim : yLim : yLim);
        grid on;
        s2.PlotBoxAspectRatio = s1.PlotBoxAspectRatio;
        
        % Extract intensity profile along horizontal centerline (y = 0)
        intensityRow = intensity(round(Ny/2), :);
        
        % Plot horizontal profile in subplot (bottom-left)
        s4 = subplot(4,3,4);
        plot(x0*1e3, intensityRow, 'k', 'LineWidth', 1.2);
        xlim([-xLim xLim]);
        xlabel('$x$ (mm)');
        ylabel('Intensity along $y = 0$');
        xticks(-xLim : xLim : xLim);
        grid on;
        s4.PlotBoxAspectRatio = s2.PlotBoxAspectRatio;
        
        %% ------- [ Detect Closest Local Minima to the Center ] -------
        % Detect approximate local minima in horizontal intensity profile

        diffRow = diff(sign(diff(intensityRow)));
        idxMinRow = find(diffRow == 2) + 1;
        
        % Detect local minima in vertical intensity profile
        diffCol = diff(sign(diff(intensityCol)));
        idxMinCol = find(diffCol == 2) + 1;
        
        % Index of the beam center
        minRowCenterIdx = round(Nx/2) + 1;
        minColCenterIdx = round(Ny/2) + 1;
        
        % Compute absolute distances from center to each local minimum
        distMinRow = abs(idxMinRow - minRowCenterIdx);
        distMinCol = abs(idxMinCol - minColCenterIdx);
        
        % Sort minima by proximity to center
        [~, sortedIdxRow] = sort(distMinRow);
        [~, sortedIdxCol] = sort(distMinCol);
        
        % Select two closest local minima to center
        twoClosestMinRow = idxMinRow(sortedIdxRow(1:2));
        twoClosestMinCol = idxMinCol(sortedIdxCol(1:2));
        
        % Retrieve intensity values at closest local minima
        intensityClosestMinRow = intensityRow(twoClosestMinRow);
        intensityClosestMinCol = intensityCol(twoClosestMinCol);
        
        %% ------- [ Combined Profile Plot with Annotated Minima ] -------

        s4 = subplot(4,3,4); hold on;
        plot(x0(twoClosestMinRow)*1000, intensityClosestMinRow, 'g*');
        
        s2 = subplot(4,3,2); hold on;
        plot(intensityClosestMinCol, y0(twoClosestMinCol)*1000, 'g*');
        
        s789 = subplot(4,3,[7 8 9]); hold on;
        
        % Horizontal profile: intensity along y = 0
        p1 = plot(x0*1e3, intensityRow, 'LineWidth', 1.2);
        p2 = plot(x0(idxMinRow)*1e3, intensityRow(idxMinRow), 'o');
        p3 = plot(x0(twoClosestMinRow)*1e3, intensityClosestMinRow, '.', 'MarkerSize', 40);
        
        % Vertical profile: intensity along x = 0
        p4 = plot(y0*1e3, intensityCol, 'LineWidth', 1.2);
        p5 = plot(y0(idxMinCol)*1e3, intensityCol(idxMinCol), 'o');
        p6 = plot(y0(twoClosestMinCol)*1e3, intensityClosestMinCol, '.', 'MarkerSize', 20);
        
        xlim([-min(xLim,yLim) min(xLim,yLim)]);
        box on;
        
        legend([p1, p2, p3, p4, p5, p6], ...
            {'Intensity along $y = 0$', ...
             'Local minima on $y = 0$', ...
             '2 closest minima to center on $y = 0$', ...
             'Intensity along $x = 0$', ...
             'Local minima on $x = 0$', ...
             '2 closest minima to center on $x = 0$'}, ...
             'Location', 'northeast');
        
        %% ------- [ Circular Mask Generation for Central Spot Isolation ] -------

        % Determine minimum intensity among closest minima
        minIntensity = min([intensityClosestMinRow, intensityClosestMinCol']);
        
        % Determine minimum radial distance from center to closest minima
        minRadius = min(abs([x0(twoClosestMinRow), y0(twoClosestMinCol)]));
        
        % Create binary circular mask
        mask = double(r < minRadius) .* double(intensity > minIntensity);
        
        % Apply mask to isolate central spot
        spotIntensity = intensity .* mask;
        
        % ------- [ Save Central Spot Matrix to File ] -------
        outputFilename = fullfile(dirData, sprintf('spotIntensity_m2=%d_m1=%d.txt', m2, m1));
        writematrix(spotIntensity, outputFilename, 'Delimiter', 'tab');
        
        %% ------- [ Visualization of Central Spot Extraction Process ] -------
        % Subplot 1: Original beam intensity

        subplot(4,3,10);
        imagesc(x0*1e+3, y0*1e+3, intensity);
        axis image;
        set(gca, 'YDir', 'normal'); colormap(hot); colorbar();
        xLim = 0.4/4; yLim = 0.3/3;
        xlim([-xLim, xLim]); ylim([-yLim, yLim]);
        xticks(-xLim:xLim:xLim); yticks(-yLim:yLim:yLim);
        xlabel('$x$ (mm)'); ylabel('$y$ (mm)');
        clim(get(gca, "CLim") .* [0 1])
        climVal = clim;
        
        % Subplot 2: Binary circular mask
        subplot(4,3,11);
        imagesc(x0*1e+3, y0*1e+3, mask);
        axis image;
        set(gca, 'YDir', 'normal'); colormap(gray); colorbar();
        xlim([-xLim, xLim]); ylim([-yLim, yLim]);
        xticks(-xLim:xLim:xLim); yticks(-yLim:yLim:yLim);
        xlabel('$x$ (mm)'); ylabel('$y$ (mm)');
        
        % Subplot 3: Extracted central spot
        subplot(4,3,12);
        imagesc(x0*1e+3, y0*1e+3, spotIntensity);
        axis image;
        set(gca, 'YDir', 'normal'); colormap(hot); colorbar();
        clim(climVal);
        xlim([-xLim, xLim]); ylim([-yLim, yLim]);
        xticks(-xLim:xLim:xLim); yticks(-yLim:yLim:yLim);
        xlabel('$x$ (mm)'); ylabel('$y$ (mm)');
        
        applyLatexToFigure(fig, 10)
        
        % Export central spot analysis figure
        fileName = sprintf('CentralSpotAnalysis_m2=%d_m1=%d.png', m2, m1);
        exportgraphics(fig, fullfile(dirImage, fileName), 'Resolution', 300);
        close(fig)
        
        % Calculate metrics
        totalPower = sum(sum(spotIntensity));
        area = sum(sum(mask));
        powerDensity = totalPower / area;
        originIntensity = intensity(round(Ny/2), round(Nx/2));
        
        analysis(m1Idx).m1 = m1;
        analysis(m1Idx).totalPower = totalPower;
        analysis(m1Idx).area = area;
        analysis(m1Idx).powerDensity = powerDensity;
        analysis(m1Idx).originIntensity = originIntensity;
        
    end
    
    %% ------- [ Relative Feature Normalization with Respect to Base m1 ] -------

    baseTotalPower = analysis(1).totalPower;
   
    basePowerDensity = analysis(1).powerDensity;
    baseOriginIntensity = analysis(1).originIntensity;
    baseArea =  analysis(1).area;
    
    relTotalPower = zeros(1, length(m1Values));
    relPowerDensity = zeros(1, length(m1Values));
    relOriginIntensity = zeros(1, length(m1Values));
    relArea = zeros(1, length(m1Values));
    
    for ii = 1:length(m1Values)
        relTotalPower(ii) = analysis(ii).totalPower / baseTotalPower;
        relPowerDensity(ii) = analysis(ii).powerDensity / basePowerDensity;
        relOriginIntensity(ii) = analysis(ii).originIntensity / baseOriginIntensity;
        relArea(ii) = analysis(ii).area / baseArea;
    end
    
    result = [[m2;m2;m2] [1*m2;3*m2;5*m2] relTotalPower' relPowerDensity' relOriginIntensity' relArea'];
    allResults = [allResults; result];

    %% ------- [ Visualization and Reporting of Normalized Beam Properties ] -------

    fig = figure('Visible', 'on');

    % Plot 1: Relative Total Power
    subplot(4,1,1);
    vals0 = relTotalPower;
    b0 = bar(m1Values, vals0, 'FaceColor', [0.9 0.5 0.2]);
    xlabel('Radial grating coefficient ($m_1$)');
    ylabel('Rel. Total Power');
    title('Normalized Total Power vs $m_1$');
    grid on;
    for ii = 1:length(m1Values)
        text(m1Values(ii), vals0(ii) + 0.01, sprintf('%.2f', vals0(ii)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    
    % Plot 2: Relative Power Density
    subplot(4,1,2);
    vals1 = relPowerDensity;
    b1 = bar(m1Values, vals1, 'FaceColor', [0.2 0.6 0.8]);
    xlabel('Radial grating coefficient ($m_1$)');
    ylabel('Rel. Power Density');
    title('Normalized Power Density vs $m_1$');
    grid on;
    for ii = 1:length(m1Values)
        text(m1Values(ii), vals1(ii) + 0.01, sprintf('%.2f', vals1(ii)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    
    % Plot 3: Relative Origin Intensity
    subplot(4,1,3);
    vals2 = relOriginIntensity;
    b2 = bar(m1Values, vals2, 'FaceColor', [0.6 0.3 0.6]);
    xlabel('Radial grating coefficient ($m_1$)');
    ylabel('Rel. Origin Intensity');
    title('Normalized Central Intensity vs $m_1$');
    grid on;
    for ii = 1:length(m1Values)
        text(m1Values(ii), vals2(ii) + 0.01, sprintf('%.2f', vals2(ii)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    
    % Plot 4: Relative Area of Central Spot
    subplot(4,1,4);
    vals3 = relArea;
    b3 = bar(m1Values, vals3, 'FaceColor', [0.3 0.7 0.3]);
    xlabel('Radial grating coefficient ($m_1$)');
    ylabel('Rel. Spot Area');
    title('Normalized Spot Area vs $m_1$');
    grid on;
    for ii = 1:length(m1Values)
        text(m1Values(ii), vals3(ii) + 0.01, sprintf('%.2f', vals3(ii)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    
    applyLatexToFigure(fig)

    % Export normalized metrics bar plot
    fileName = sprintf('NormalizedMetrics_m2=%d.png', m2);
    exportgraphics(gcf, fullfile(dirImage, fileName), 'Resolution', 300);
    close(fig)

    % Console Output
    disp('Comparison results:');
    fprintf('%5s | %12s | %12s | %12s | %12s\n', 'm1', 'Rel. Total Power', 'Rel. Power', 'Rel. Intensity', 'Rel. Area');
    fprintf('%s\n', repmat('-', 1, 68));
    for ii = 1:length(m1Values)
        fprintf(' %3d  |     %7.4f   |     %7.4f   |     %7.4f   |     %7.4f\n', ...
            m1Values(ii), ...
            relTotalPower(ii), ...
            relPowerDensity(ii), ...
            relOriginIntensity(ii), ...
            relArea(ii));
    end
    
    %% ------- [ Export of Normalized Results to Text File ] -------
    % Export text report

    outputFilename = fullfile(dirResult, sprintf('spot_metrics_summary_sim_m2=%d.txt', m2));
    
    fid = fopen(outputFilename, 'w');
    fprintf(fid, 'Comparison of Normalized Quantities for Different m1 (m2 = %d)\n', m2);
    fprintf(fid, '%5s | %15s | %18s | %20s | %15s\n', ...
            'm1', 'Rel. Total Power', 'Rel. Power Density', ...
            'Rel. Origin Intensity', 'Rel. Spot Area');
    fprintf(fid, '%s\n', repmat('-', 1, 85));
    
    for ii = 1:length(m1Values)
        fprintf(fid, '%5d | %15.4f | %18.4f | %20.4f | %15.4f\n', ...
            m1Values(ii), ...
            relTotalPower(ii), ...
            relPowerDensity(ii), ...
            relOriginIntensity(ii), ...
            relArea(ii));
    end
    
    fclose(fid);
    fprintf('Results saved to: %s\n', outputFilename);
    close all
end

% Save aggregate results table
T = array2table(allResults, ...
    'VariableNames', {'m2', 'm1', 'RelTotalPower', 'RelPowerDensity', 'RelOriginIntensity', 'RelArea'});
writetable(T, fullfile(dirResult, 'spot_metrics_summary_sim.txt'), 'Delimiter', '\t');

%% ------- [ Load allresults if Empty ] -------
% If 'allresults' is empty or undefined, load it from saved table file

if ~exist('allresults', 'var') || isempty(allResults)
    fprintf('Loading allresults from saved table...\n');
    T = readtable(fullfile(dirResult, 'spot_metrics_summary_sim.txt'), 'Delimiter', '\t');
    allResults = table2array(T);
end

%% ------- [ Grouped Bar Plot of Relative Power Densities ] -------
% Visualize relative power density for m1 = [1×m2, 3×m2, 5×m2] across m2 values

target_m2 = [5 10];
filteredRows = ismember(allResults(:,1), target_m2);
filteredResults = allResults(filteredRows, :);
m2Unique = unique(filteredResults(:,1));
m2Num = length(m2Unique);

% Initialize matrix for relative power density values
relPowerDensity = zeros(m2Num, 3);

for ii = 1:m2Num
    idx = allResults(:,1) == m2Unique(ii);
    groupData = allResults(idx, :);
    [~, sortIdx] = sort(groupData(:,2));
    sortedGroup = groupData(sortIdx, :);
    relPowerDensity(ii, :) = sortedGroup(:,4)';
end

% Create grouped bar chart
figure('Position', [200 406 1022 383]);
b = bar(relPowerDensity', 'grouped');

colors = [
    0.1216 0.4667 0.7059
    1.0000 0.4980 0.0549
    0.1725 0.6275 0.1725
    0.8392 0.1529 0.1569
    0.5804 0.4039 0.7412
    0.5490 0.3373 0.2941
    0.8902 0.4667 0.7608
    0.4980 0.4980 0.4980
    0.7373 0.7412 0.1333
    0.0902 0.7451 0.8118
    1.0000 0.8431 0.0000
];

for k = 1:m2Num
    b(k).FaceColor = colors(k,:);
end

xlabel('$m_1 = \left[1\times m_2, 3\times m_2, 5\times m_2\right]$');
ylabel('Relative Power Density');
title('Grouped Bar Chart of Relative Power Density for Each $m_2$');
xticklabels({'$1\times m_2$','$3\times m_2$','$5\times m_2$'});
legend(arrayfun(@(m) sprintf('$m_2 = %d$', m), m2Unique, 'UniformOutput', false), ...
    'Location', 'northeastoutside');
grid on;
yticks(0:0.05:1);

applyLatexToFigure(gcf)

% Export grouped bar charts
exportgraphics(gcf, fullfile(dirImage, 'rel_power_density.pdf'));
exportgraphics(gcf, fullfile(dirImage, 'rel_power_density.png'), 'Resolution', 300);
close(gcf)

%% ------- [ Grouped Bar Plot of Total Power ] -------
% Visualize relative total power for m1 = [1×m2, 3×m2, 5×m2] across m2 values

relTotalPower = zeros(m2Num, 3);

for ii = 1:m2Num
    idx = allResults(:,1) == m2Unique(ii);
    groupData = allResults(idx, :);
    [~, sortIdx] = sort(groupData(:,2));
    sortedGroup = groupData(sortIdx, :);
    relTotalPower(ii, :) = sortedGroup(:,3)';
end

figure('Position', [200 406 1022 383]);
b = bar(relTotalPower', 'grouped');

for k = 1:m2Num
    b(k).FaceColor = colors(k,:);
end

xlabel('$m_1 = \left[1\times m_2, 3\times m_2, 5\times m_2\right]$');
ylabel('Relative Total Power');
title('Grouped Bar Chart of Relative Total Power for Each $m_2$');
xticklabels({'$1\times m_2$','$3\times m_2$','$5\times m_2$'});
legend(arrayfun(@(m) sprintf('$m_2 = %d$', m), m2Unique, 'UniformOutput', false), ...
    'Location', 'northeastoutside');
grid on;
yticks(0:0.05:1);

applyLatexToFigure(gcf)

exportgraphics(gcf, fullfile(dirImage, 'rel_total_power.pdf'));
exportgraphics(gcf, fullfile(dirImage, 'rel_total_power.png'), 'Resolution', 300);
close(gcf)

%% ------- [ Grouped Bar Plot of Origin Intensity ] -------
% Visualize relative origin intensity for m1 = [1×m2, 3×m2, 5×m2] across m2 values

relOriginIntensity = zeros(m2Num, 3);

for ii = 1:m2Num
    idx = allResults(:,1) == m2Unique(ii);
    groupData = allResults(idx, :);
    [~, sortIdx] = sort(groupData(:,2));
    sortedGroup = groupData(sortIdx, :);
    relOriginIntensity(ii, :) = sortedGroup(:,5)';
end

figure('Position', [200 406 1022 383]);
b = bar(relOriginIntensity', 'grouped');

for k = 1:m2Num
    b(k).FaceColor = colors(k,:);
end

xlabel('$m_1 = \left[1\times m_2, 3\times m_2, 5\times m_2\right]$');
ylabel('Relative Origin Intensity');
title('Grouped Bar Chart of Relative Origin Intensity for Each $m_2$');
xticklabels({'$1\times m_2$','$3\times m_2$','$5\times m_2$'});
legend(arrayfun(@(m) sprintf('$m_2 = %d$', m), m2Unique, 'UniformOutput', false), ...
    'Location', 'northeastoutside');
grid on;
yticks(0:0.05:1);

applyLatexToFigure(gcf)

exportgraphics(gcf, fullfile(dirImage, 'rel_origin_intensity.pdf'));
exportgraphics(gcf, fullfile(dirImage, 'rel_origin_intensity.png'), 'Resolution', 300);
close(gcf)

%% ------- [ Grouped Bar Plot of Spot Area ] -------
% Visualize relative spot area for m1 = [1×m2, 3×m2, 5×m2] across m2 values

relArea = zeros(m2Num, 3);

for ii = 1:m2Num
    idx = allResults(:,1) == m2Unique(ii);
    groupData = allResults(idx, :);
    [~, sortIdx] = sort(groupData(:,2));
    sortedGroup = groupData(sortIdx, :);
    relArea(ii, :) = sortedGroup(:,6)';
end

figure('Position', [200 406 1022 383]);
b = bar(relArea', 'grouped');

for k = 1:m2Num
    b(k).FaceColor = colors(k,:);
end

xlabel('$m_1 = \left[1\times m_2, 3\times m_2, 5\times m_2\right]$');
ylabel('Relative Spot Area');
title('Grouped Bar Chart of Relative Spot Area for Each $m_2$');
xticklabels({'$1\times m_2$','$3\times m_2$','$5\times m_2$'});
legend(arrayfun(@(m) sprintf('$m_2 = %d$', m), m2Unique, 'UniformOutput', false), ...
    'Location', 'northeastoutside');
grid on;
yticks(0:0.05:1);

applyLatexToFigure(gcf)

exportgraphics(gcf, fullfile(dirImage, 'rel_area.pdf'));
exportgraphics(gcf, fullfile(dirImage, 'rel_area.png'), 'Resolution', 300);
close(gcf)

%% ------- [ Line Plots of Normalized Quantities vs m1/m2 ] -------
% Plot four normalized quantities as line plots with markers

ratios = [1, 3, 5];
figure('Position', [100, 100, 1000, 700]);

% Plot 1: Relative Total Power
subplot(2,2,1);
plot(ratios, relTotalPower', 'o', 'LineWidth', 1.5);
xlabel('$m_1$/$m_2$');
ylabel('Rel. Total Power');
title('Normalized Total Power');
grid on;
legend(arrayfun(@(m) sprintf('$m_2 = %d$', m), m2Unique, 'UniformOutput', false), ...
    'Location', 'best');

% Plot 2: Relative Power Density
subplot(2,2,2);
plot(ratios, relPowerDensity', 'o', 'LineWidth', 1.5);
xlabel('$m_1$/$m_2$');
ylabel('Rel. Power Density');
title('Normalized Power Density');
grid on;
legend(arrayfun(@(m) sprintf('$m_2 = %d$', m), m2Unique, 'UniformOutput', false), ...
    'Location', 'best');

% Plot 3: Relative Origin Intensity
subplot(2,2,3);
plot(ratios, relOriginIntensity', 'o', 'LineWidth', 1.5);
xlabel('$m_1$/$m_2$');
ylabel('Rel. Origin Intensity');
title('Normalized Central Intensity');
grid on;
legend(arrayfun(@(m) sprintf('$m_2 = %d$', m), m2Unique, 'UniformOutput', false), ...
    'Location', 'best');

% Plot 4: Relative Area
subplot(2,2,4);
plot(ratios, relArea', 'o', 'LineWidth', 1.5);
xlabel('$m_1$/$m_2$');
ylabel('Rel. Spot Area');
title('Normalized Spot Area');
grid on;
legend(arrayfun(@(m) sprintf('$m_2 = %d$', m), m2Unique, 'UniformOutput', false), ...
    'Location', 'best');

applyLatexToFigure(gcf)

% Export line plots
exportgraphics(gcf, fullfile(dirImage, 'lineplots_all_quantities.pdf'), 'ContentType', 'vector');
close(gcf)

%% ------- [ Process All Lab Images ] -------
% Find radial symmetry center and extract spot metrics
% for all images named m1=#.png in ../labImage

% Get list of all images in the folder
imageFolder = fullfile('..', 'labImage');
imageFiles = dir(fullfile(imageFolder, 'm1=*.png'));

if isempty(imageFiles)
    warning('No image files found in: %s', imageFolder);
    return;
end

fprintf('Found %d image(s) in %s\n', length(imageFiles), imageFolder);

% Initialize table data
m1Values = [];
spotPowerDensity = [];
spotMaxIntensity = [];

% Process each image
for fileIdx = 1:length(imageFiles)
    fileName = imageFiles(fileIdx).name;
    filePath = fullfile(imageFolder, fileName);
    
    % Extract m1 from filename
    tokens = regexp(fileName, 'm1=(\d+)', 'tokens');
    if isempty(tokens)
        warning('Could not extract m1 from filename: %s', fileName);
        continue;
    end
    m1 = str2double(tokens{1}{1});
    
    fprintf('\nProcessing m1 = %d...\n', m1);
    
    % Read image
    originalImage = imread(filePath);
    fprintf('  Image loaded: %s (%d x %d x %d)\n', fileName, size(originalImage,1), size(originalImage,2), size(originalImage,3));
    
    % Convert to grayscale
    if size(originalImage, 3) == 3
        grayImage = double(rgb2gray(originalImage));
    else
        grayImage = double(originalImage);
    end
    grayImage = grayImage / max(grayImage(:));
    
    % Find center using grayscale image
    [xCenter, yCenter] = findRadialSymmetryCenter(grayImage);
    fprintf('  Center from grayscale: (%.2f, %.2f)\n', xCenter, yCenter);
    
    % Extract spot features from grayscale
    [hasSpot, estRadius, spotArea, spotSum, spotPowerDensityVal, spotMaxIntensityVal] = ...
        spotDetection(grayImage, xCenter, yCenter);
    
    % Store values for table
    m1Values = [m1Values; m1];
    spotPowerDensity = [spotPowerDensity; spotPowerDensityVal];
    spotMaxIntensity = [spotMaxIntensity; spotMaxIntensityVal];
    
    % Display results
    fig = figure('Color', 'white');
    
    imshow(originalImage);
    hold on;
    % Plot center from green channel
    plot(xCenter, yCenter, 'r+', 'MarkerSize', 20, 'LineWidth', 3);
    plot(xCenter, yCenter, 'ro', 'MarkerSize', 25, 'LineWidth', 2);
    % Add cross-hair lines
    xlim([0.5, size(originalImage, 2) + 0.5]);
    ylim([0.5, size(originalImage, 1) + 0.5]);
    plot([1, size(originalImage, 2)], [yCenter, yCenter], 'r--', 'LineWidth', 0.5);
    plot([xCenter, xCenter], [1, size(originalImage, 1)], 'r--', 'LineWidth', 0.5);
    
    % Draw spot radius circle
    if hasSpot
        theta = linspace(0, 2*pi, 100);
        xCircle = xCenter + estRadius * cos(theta);
        yCircle = yCenter + estRadius * sin(theta);
        plot(xCircle, yCircle, 'g-', 'LineWidth', 2);
    end
    hold off;

    title(sprintf('Center: (%.1f, %.1f)', xCenter, yCenter), ...
          'FontSize', 12, 'FontName', 'Times New Roman');
    
    % Export figure
    exportgraphics(fig, fullfile(dirImage, sprintf('center_detection_m1=%d.png', m1)), 'Resolution', 300);
    close(fig);
end

%% ------- [ Normalize to m1 = 10 ] -------
% Find the reference value for m1 = 10
idxRef = find(m1Values == 10);
if isempty(idxRef)
    warning('m1 = 10 not found. Using first value as reference.');
    refPowerDensity = spotPowerDensity(1);
    refMaxIntensity = spotMaxIntensity(1);
else
    refPowerDensity = spotPowerDensity(idxRef);
    refMaxIntensity = spotMaxIntensity(idxRef);
    fprintf('\nReference (m1 = 10):\n');
    fprintf('  SpotPowerDensity: %.4f\n', refPowerDensity);
    fprintf('  SpotMaxIntensity: %.4f\n', refMaxIntensity);
end

% Normalize all values
spotPowerDensityNorm = spotPowerDensity / refPowerDensity;
spotMaxIntensityNorm = spotMaxIntensity / refMaxIntensity;

fprintf('\nNormalized values:\n');
for i = 1:length(m1Values)
    fprintf('  m1 = %d: PowerDensity = %.4f, MaxIntensity = %.4f\n', ...
            m1Values(i), spotPowerDensityNorm(i), spotMaxIntensityNorm(i));
end

%% ------- [ Create Results Table (Normalized Only) ] -------
% Create table with only normalized values
resultsTable = table(m1Values, spotPowerDensityNorm, spotMaxIntensityNorm, ...
                     'VariableNames', {'m1', 'SpotPowerDensity', 'SpotMaxIntensity'});

% Display the table
disp(resultsTable);

writetable(resultsTable, fullfile(dirResult, 'spot_metrics_summary_lab.txt'), 'Delimiter', '\t');
fprintf('\nTable saved to: %s\n', fullfile(dirResult, 'spot_metrics_summary_lab.txt'));

%% ------- [ Print Summary Table ] -------
% Display results in a formatted table

fprintf('\n');
fprintf('================================================================================\n');
fprintf('                   NORMALIZED SPOT ANALYSIS SUMMARY TABLE\n');
fprintf('================================================================================\n');
fprintf('  m1    SpotPowerDensity    SpotMaxIntensity\n');
fprintf('--------------------------------------------------------------------------------\n');

for i = 1:length(m1Values)
    fprintf('  %2d         %8.4f            %8.4f\n', ...
            m1Values(i), ...
            spotPowerDensityNorm(i), ...
            spotMaxIntensityNorm(i));
end

fprintf('================================================================================\n');
fprintf('Normalized to m1 = 10 (reference)\n');

%% ------- [ Compare Sim and Lab Power Density ] -------
% Read sim and lab tables and plot bar chart for m1 = 10, 25, 50

% Read sim table
simFile = fullfile(dirResult, 'spot_metrics_summary_sim.txt');
if exist(simFile, 'file')
    simTable = readtable(simFile, 'Delimiter', '\t');
    fprintf('Sim table loaded: %s\n', simFile);
else
    error('Sim table not found: %s', simFile);
end

% Read lab table
labFile = fullfile(dirResult, 'spot_metrics_summary_lab.txt');
if exist(labFile, 'file')
    labTable = readtable(labFile, 'Delimiter', '\t');
    fprintf('Lab table loaded: %s\n', labFile);
else
    error('Lab table not found: %s', labFile);
end

% Display tables
disp('Sim Table:');
disp(simTable);
disp('Lab Table:');
disp(labTable);

% Extract data for specific m1 values (m1 = 10, 25, 50)
targetM1 = [10, 30, 50];

% Initialize arrays
simPowerDensity = zeros(length(targetM1), 1);
labPowerDensity = zeros(length(targetM1), 1);

for i = 1:length(targetM1)
    m1 = targetM1(i);
    
    % Find in sim table (m1 column)
    idxSim = find(simTable.m1 == m1);
    if ~isempty(idxSim)
        simPowerDensity(i) = simTable.RelPowerDensity(idxSim);
    else
        warning('m1 = %d not found in sim table', m1);
        simPowerDensity(i) = NaN;
    end
    
    % Find in lab table (m1 column)
    idxLab = find(labTable.m1 == m1);
    if ~isempty(idxLab)
        labPowerDensity(i) = labTable.SpotPowerDensity(idxLab);
    else
        warning('m1 = %d not found in lab table', m1);
        labPowerDensity(i) = NaN;
    end
end

% Create figure with two subplots
fig = figure('Color', 'white', 'Position', [123.00 454.00 850.00 200.00]);

% ----- Subplot 1: Simulation -----
subplot(1, 2, 1);
b1 = bar(targetM1, simPowerDensity, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k', 'LineWidth', 0.5);
xlabel('$m_1$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Relative Power Density', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Simulation', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
grid minor;
ylim([0, 1]);
yticks(0:0.2:1);
xlim([-5 65])

% Add value labels
for i = 1:length(targetM1)
    text(targetM1(i), simPowerDensity(i) + 0.03, sprintf('%.2f', simPowerDensity(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
end

% ----- Subplot 2: Experiment -----
subplot(1, 2, 2);
b2 = bar(targetM1, labPowerDensity, 'FaceColor', [0.9 0.5 0.2], 'EdgeColor', 'k', 'LineWidth', 0.5);
xlabel('$m_1$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Relative Power Density', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Experiment', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
grid minor;
ylim([0, 1]);
yticks(0:0.2:1);
xlim([-5 65])

% Add value labels
for i = 1:length(targetM1)
    text(targetM1(i), labPowerDensity(i) + 0.03, sprintf('%.2f', labPowerDensity(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
end

% Apply LaTeX formatting
applyLatexToFigure(fig, 12);

% Export figure
exportgraphics(fig, fullfile(dirImage, 'sim_vs_lab_power_density.pdf'));
exportgraphics(fig, fullfile(dirImage, 'sim_vs_lab_power_density.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', fullfile(dirImage, 'sim_vs_lab_power_density.pdf'));
close(fig);

%% ------- [ Local Functions ] -------

function [xCenter, yCenter] = findRadialSymmetryCenter(image)
% findRadialSymmetryCenter - Estimates the center of 180° radial symmetry in a grayscale image.
%
% Syntax:
%   [xCenter, yCenter] = findRadialSymmetryCenter(grayDoubleImageReduced)
%
% Description:
%   This function searches for the point in the image around which the 180-degree
%   rotated version of the image most closely matches the original. It first performs
%   a coarse grid search using integer-pixel shifts and then refines the estimated
%   center using subpixel translation to achieve higher accuracy.
%
%   The method assumes that the image contains a radially symmetric pattern centered
%   near the image center (e.g., Poisson-Arago spot), and returns the coordinates of
%   the estimated symmetry center with subpixel precision.
%
% Input:
%   grayDoubleImageReduced - A 2D grayscale image of type double, preprocessed and normalized.
%
% Outputs:
%   xCenter - Estimated X-coordinate of the radial symmetry center (may be fractional)
%   yCenter - Estimated Y-coordinate of the radial symmetry center (may be fractional)

    [rows, cols] = size(image);
    delta = 10;
    cx = round(cols / 2);
    cy = round(rows / 2);
    minCost = inf;
    bestX = cx;
    bestY = cy;

    margin = 5;
    mask = ones(rows, cols);
    mask(1:margin, :) = 0;
    mask(end-margin:end, :) = 0;
    mask(:, 1:margin) = 0;
    mask(:, end-margin:end) = 0;
    
    rotated = imrotate(image, 180, 'bilinear', 'crop');

    for dx = -delta:delta
        for dy = -delta:delta
            ox = cx + dx;
            oy = cy + dy;
            shifted = circshift(rotated, [2*dy, 2*dx]);
            diff = abs(image - shifted) .* mask;
            score = sum(diff(:));
            if score < minCost
                minCost = score;
                bestX = ox;
                bestY = oy;
            end
        end
    end

    fineStep = 0.1;
    fineRange = -0.5:fineStep:0.5;
    minFineCost = inf;
    refinedX = bestX;
    refinedY = bestY;

    for fx = fineRange
        for fy = fineRange
            shiftX = 2 * (bestX - cx + fx);
            shiftY = 2 * (bestY - cy + fy);
            translated = imtranslate(rotated, [shiftX, shiftY], 'OutputView', 'same', 'Method', 'linear');
            diffFine = abs(image - translated) .* mask;
            fineScore = sum(diffFine(:));
            if fineScore < minFineCost
                minFineCost = fineScore;
                refinedX = bestX + fx;
                refinedY = bestY + fy;
            end
        end
    end

    xCenter = refinedX;
    yCenter = refinedY;
end

function [hasSpot, estRadius, spotArea, spotSum, spotPowerDensity, spotMaxIntensity] = spotDetection(image, xCenter, yCenter)
% spotDetection - Detects a central optical spot and computes its properties.
%
% Syntax:
%   [hasSpot, estRadius, spotArea, spotSum, spotPower] = ...
%       spotDetection(image, xCenter, yCenter)
%
% Description:
%   This function analyzes the radial intensity distribution around a given
%   center point (xCenter, yCenter) in a grayscale image to detect the
%   presence of a central bright spot (e.g., Poisson-Arago spot). If a spot
%   is detected, it estimates its radius, total intensity, and average power.
%
% Inputs:
%   image    - 2D grayscale image (type double), preferably normalized [0,1]
%   xCenter  - X-coordinate of the candidate spot center
%   yCenter  - Y-coordinate of the candidate spot center
%
% Outputs:
%   hasSpot    - Logical flag (true/false) indicating if a central spot is detected
%   estRadius  - Estimated radius of the bright spot (in pixels)
%   spotArea   - Number of pixels within the estimated spot region
%   spotSum    - Total intensity (sum) inside the spot
%   spotPower  - Mean intensity (average power) of the spot

    [rows, cols] = size(image);
    [X, Y] = meshgrid(1:cols, 1:rows);
    R = sqrt((X - xCenter).^2 + (Y - yCenter).^2);
    
    % ----- Build initial radial intensity profile -----
    maxR = min([rows, cols]) / 25;
    radii = 1:round(maxR);
    radialProfile = zeros(size(radii));
    
    for idx = 1:length(radii)
        mask = abs(R - radii(idx)) < 0.5;
        radialProfile(idx) = mean(image(mask));
    end
    
    profileLen = length(radialProfile);
    
    % --- Compute outerVal dynamically and safely ---
    outerIdx = max(round(profileLen * 0.6), profileLen - 3);
    outerVal = mean(radialProfile(outerIdx:end));
    
    % --- Compute innerVal based on brightness level ---
    if outerVal < 0.5
        % Use first ~30% of profile, or at least 1 point, but avoid overlap with outer region
        innerIdx = min(round(profileLen * 0.3), outerIdx - 1);
        innerIdx = max(innerIdx, 1);  % Ensure at least one point is used
        innerVal = mean(radialProfile(1:innerIdx));
    else
        innerVal = radialProfile(1);  % Use only the center when brightness is sufficient
    end
    
    dropRatio = innerVal / outerVal;
    threshold = 1.05 + 0.1 * (outerVal < 0.5);
    
    if (dropRatio < threshold)
        hasSpot = false;
        estRadius = NaN;
        spotArea = 0;
        spotSum = 0;
        spotPowerDensity = 0;
        spotMaxIntensity = 0;
        return;
    end
    
    % ----- Smooth and analyze radial profile -----
    smoothWin = min(profileLen, max(3, round(profileLen * 0.1)));
    smoothedProfile = movmean(radialProfile, smoothWin);
    normProfile = smoothedProfile / max(smoothedProfile);
    gradProfile = gradient(normProfile);
    
    % ----- Estimate spot radius -----
    gradientThreshold = -0.05;
    dropIndex = find(gradProfile < gradientThreshold, 1);
    estRadius = min(dropIndex, 15);
    hasSpot = true;
    
    % ----- Compute spot metrics -----
    spotMask = R <= estRadius;
    spotArea = nnz(spotMask);
    spotSum = sum(image(spotMask));
    spotPowerDensity = spotSum / spotArea;
    spotMaxIntensity = max(image(spotMask));    
end

function applyLatexToFigure(figHandle, fontSize)
    % Apply LaTeX interpreter and Times New Roman font to all text elements in a figure
    %
    % Inputs:
    %   figHandle - (optional) Figure handle (default: gcf)
    %   fontSize  - (optional) Font size for text and axes (default: 12)
    %
    % Usage:
    %   applyLatexToFigure();           % Apply to current figure with font size 12
    %   applyLatexToFigure(gcf);        % Apply to current figure with font size 12
    %   applyLatexToFigure(fig, 14);    % Apply to specific figure with font size 14
    
    % Set default figure handle
    if nargin < 1
        figHandle = gcf;
    end
    
    % Set default font size
    if nargin < 2
        fontSize = 12;
    end
    
    % Ensure figure exists
    if ~isvalid(figHandle)
        error('Invalid figure handle.');
    end
    
    % Activate the figure
    figure(figHandle);
    
    % ----- Apply to Axes and All Children -----
    % Get all axes in the figure (including subplots)
    allAxes = findall(figHandle, 'Type', 'axes');
    
    for i = 1:length(allAxes)
        ax = allAxes(i);
        
        % Set axes properties
        set(ax, 'FontName', 'Times New Roman', 'FontSize', fontSize);
        set(ax, 'TickLabelInterpreter', 'latex');
        set(ax, 'DefaultTextInterpreter', 'latex');
        
        % ----- Axis Labels -----
        xlabel(ax, get(get(ax, 'XLabel'), 'String'), ...
               'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize);
        ylabel(ax, get(get(ax, 'YLabel'), 'String'), ...
               'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize);
        zlabel(ax, get(get(ax, 'ZLabel'), 'String'), ...
               'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize);
        
        % ----- Title -----
        title(ax, get(get(ax, 'Title'), 'String'), ...
              'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize + 2);
        
        % ----- Legend (if exists) -----
        leg = get(ax, 'Legend');
        if ~isempty(leg) && isvalid(leg)
            set(leg, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize - 2);
        end
        
        % ----- Colorbar (if exists) -----
        cb = get(ax, 'Colorbar');
        if ~isempty(cb) && isvalid(cb)
            try
                set(cb, 'FontName', 'Times New Roman', 'FontSize', fontSize - 2);
                set(cb, 'TickLabelInterpreter', 'latex');
            catch
                % If colorbar doesn't support these properties, skip
            end
            
            % Colorbar label
            try
                cbLabel = get(cb, 'Label');
                if ~isempty(cbLabel) && isvalid(cbLabel)
                    set(cbLabel, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize);
                end
            catch
                % If colorbar label doesn't exist or can't be modified, skip
            end
        end
    end
    
    % ----- Text Annotations (if any) -----
    allText = findall(figHandle, 'Type', 'text');
    for i = 1:length(allText)
        try
            set(allText(i), 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize);
        catch
            % If text doesn't support these properties, skip
        end
    end
    
    % ----- Suptitle (if exists) -----
    allSuptitles = findall(figHandle, 'Type', 'text', 'Tag', 'suptitle');
    for i = 1:length(allSuptitles)
        try
            set(allSuptitles(i), 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize + 4);
        catch
            % If suptitle doesn't support these properties, skip
        end
    end
    
    % ----- Handle sgtitle (Figure title created by sgtitle) -----
    allSgtitles = findall(figHandle, 'Type', 'text', 'Tag', 'sgtitle');
    for i = 1:length(allSgtitles)
        try
            set(allSgtitles(i), 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize + 4);
        catch
            % If sgtitle doesn't support these properties, skip
        end
    end
    
    % ----- Figure Title (Name) -----
    set(figHandle, 'Name', get(figHandle, 'Name'), 'NumberTitle', 'off');
end
