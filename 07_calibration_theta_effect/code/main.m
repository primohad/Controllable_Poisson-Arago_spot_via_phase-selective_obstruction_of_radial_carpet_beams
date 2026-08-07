% ----------------------------------------------------------
% Simulate effect of angular misalignment between binary
% amplitude grating (m2) and hybrid grating (m1).
% Evaluates Poisson-Arago on-axis intensity at detection plane.
%
% Workflow: Build SLM field -> Propagate to second grating ->
%           Sweep rotation angle -> Apply amplitude grating ->
%           Propagate through lens to detection plane -> Record intensity
%
% Output: ../data/ (per-m1 tables), ../image/ (per-m1 plots),
%         ../movie/ (propagation videos)
%
% Notes: m1 = [5,10,30,50], m2 = 10 fixed.
%        Detection plane at z = f + 1 cm.
% ----------------------------------------------------------

%% ------- [ Environment Setup ] -------
% Clear workspace, close figures, compact display

clc;
clear; close all; format compact; format bank;

%% ------- [ Working Directory Sync ] -------
% Auto-change current folder to the active script’s location.

disp('Auto change folder enabled.');
activeFile = matlab.desktop.editor.getActiveFilename;
if ~isempty(activeFile)
    cd(fileparts(activeFile));
    fprintf('Changed folder to: %s\n', pwd);
end

%% ------- [ Output Paths ] -------
% Define and create output directories if they do not exist

parentDir = fullfile(pwd, '..');

dirData  = fullfile(parentDir, 'data');
dirImage = fullfile(parentDir, 'image');
dirMovie = fullfile(parentDir, 'movie');

% Create directories safely
dirPaths = {dirData, dirImage, dirMovie};
for k = 1:length(dirPaths)
    if ~exist(dirPaths{k}, 'dir')
        mkdir(dirPaths{k});
    end
end


%% ------- [ Physical Constants ] -------
% Define complex unit, wavelength, and wavenumber

i = complex(0, 1);
wl = 532e-9;    % Wavelength in meters (532 nm)
k = 2*pi/wl;    % Wavenumber [1/m]

%% ------- [ Spatial Coordinates ] -------
% Build Cartesian grid and corresponding polar coordinates

Nx = 2^12 + 1;      % Grid size along x
Ny = 2^12 + 1;      % Grid size along y

lx = 1.5e-2;        % Half-window size along x [m]
ly = 1.5e-2;        % Half-window size along y [m]

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

wx = 1.422e-2;      % SLM width [m]
wy = 1.066e-2;      % SLM height [m]

uPlane = zeros(size(x));
maskSLM = (x > -wx/2 & x < wx/2) & (y > -wy/2 & y < wy/2);
uPlane(maskSLM) = 1;  % Unit amplitude over SLM aperture

%% ------- [ Grating Parameters ] -------
% Binary radial structure with amplitude and phase modulation

Va = 0.1;       % Amplitude visibility (0 < Va < 1)
Vp = pi/2;      % Phase modulation amplitude [rad]
m2 = 10;        % Number of spokes (second grating, fixed)

%% ------- [ First Pass: Global Max Intensity ] -------
% Propagate without visualization to find global maximum for normalization
fprintf('First pass: finding global maximum intensity...\n');

globalMaxIntensity = 0;
zToSecondGrating = 2.3;    % Distance to second grating [m]
f = 20e-2;                 % Focal length of lens [m]
zToLens = 5e-2;            % Distance from second grating to lens [m]
zToDetectionPlane = f + 1e-2;  % Detection plane after lens [m]

for m1 = [5, 10, 30, 50]
    fprintf('  m1 = %d (first pass)\n', m1);
    
    % Create first hybrid radial grating
    t1 = (0.5 * (1 + Va * sign(cos(m1 * tet)))) .* ...
         exp(1i * Vp * sign(cos(m1 * tet)));
    uGratingOnSLM = t1 .* uPlane;
    
    % Propagate from SLM to second grating plane
    uToSecondGrating = ifft2(fft2(uGratingOnSLM) .* ...
        exp(-1i * pi * wl * zToSecondGrating * (fx.^2 + fy.^2)));
    uToSecondGrating = uToSecondGrating ./ sqrt(max(max(abs(uToSecondGrating).^2)));

    % Angular deviation sweep (0 to pi/m2 in 0.25 degree steps)
    angularDeviationRad = deg2rad(0:0.25:rad2deg(pi/m2));
    
    for ii = 1:length(angularDeviationRad)
        % Apply binary amplitude grating with angular shift
        t2 = 0.5 * (1 + sign(cos(m2 * (tet + angularDeviationRad(ii)))));
        uAtSecondGrating = uToSecondGrating .* t2;
        
        % Propagate to lens and apply thin lens phase
        lensPhase = exp(-1i * k * r.^2 / (2 * f));
        uAtLens = lensPhase .* ifft2(fft2(uAtSecondGrating) .* ...
            exp(-1i * pi * wl * zToLens * (fx.^2 + fy.^2)));
        
        % Propagate to detection plane
        uAtDetectionPlane = ifft2(fft2(uAtLens) .* ...
            exp(-1i * pi * wl * zToDetectionPlane * (fx.^2 + fy.^2)));
        
        % Update global maximum if current plane has higher intensity
        currentMax = max(abs(uAtDetectionPlane(:)).^2);
        if currentMax > globalMaxIntensity
            globalMaxIntensity = currentMax;
        end
    end
end

cbarMax = ceil(globalMaxIntensity);

fprintf('Global maximum intensity: %.4f\n', globalMaxIntensity);
fprintf('Ceiling value for colorbar: %d\n', ceil(globalMaxIntensity));

%% ------- [ Second Pass: Full Simulation ] -------
% Run full simulation with fixed colorbar limits
fprintf('Second pass: running full simulation...\n');

for m1 = [5, 10, 30, 50]
    fprintf('Processing m1 = %d\n', m1);

    % Create first hybrid radial grating
    t1 = (0.5 * (1 + Va * sign(cos(m1 * tet)))) .* ...
         exp(1i * Vp * sign(cos(m1 * tet)));
    
    % Field on SLM after applying the first grating
    uGratingOnSLM = t1 .* uPlane;

    % ------- [ Propagation to Second Grating | Angular Spectrum ] -------
    % Free-space propagation over zToSecondGrating using angular spectrum method.
    zToSecondGrating = 230e-2;
    uToSecondGrating = ifft2(fft2(uGratingOnSLM).*exp(-i.*pi.*wl*zToSecondGrating.*(fx.^2+fy.^2)));
    uToSecondGrating = uToSecondGrating ./ sqrt(max(max(abs(uToSecondGrating).^2)));

    % Angular deviation (theta) between the amplitude grating (m2) and the phase grating (m1)
    angularDeviationDeg = 0:0.25:rad2deg(pi/m2);
    angularDeviationRad = deg2rad(angularDeviationDeg);

    % Create figure with three panels
    fig = figure('Name', sprintf('m1 = %d', m1), 'NumberTitle', 'off', ...
                 'Units', 'normalized', 'Position', [0.07, 0.11, 0.88, 0.77], ...
                 'Color', 'white');
    s1 = subplot(2, 2, 1);  % Phase at second grating
    s2 = subplot(2, 2, 3);  % Intensity at second grating
    s3 = subplot(2, 2, [2, 4]);  % Intensity at detection plane
    
    onAxisIntensity = zeros(size(angularDeviationDeg));
    maxIntensity = zeros(size(angularDeviationDeg));
    
    % Initialize video writer for saving propagation frames
    videoFilename = sprintf('propagation_m1=%d.avi', m1);
    videoPath = fullfile(dirMovie, videoFilename);
    videoWriter = VideoWriter(videoPath);
    videoWriter.FrameRate = 5;
    open(videoWriter);

    for ii = 1:length(angularDeviationRad)
        % Apply binary amplitude radial grating with angular shift
        t2 = 0.5 * (1 + sign(cos(m2 * (tet + angularDeviationRad(ii)))));
        
        uAtSecondGrating = uToSecondGrating .* t2;
        
        % Show intensity pattern at the second grating plane
        subplot(s1);
        imagesc(x0*1000,y0*1000,abs(uAtSecondGrating).^2)
        colormap(s1,'hot'); axis image
        xlim(round(1000*[-wx/2,wx/2])); xticks(round(1000*[-wx/2,0,wx/2]))
        ylim(round(1000*[-wy/2,wy/2])); yticks(round(1000*[-wy/2,0,wy/2]))
        set(s1,'YDir','normal')
        colorbar(s1,'Ticks',[0,ceil(max(max(abs(uAtSecondGrating).^2)))], ...
            'TickLabels',{'0',ceil(max(max(abs(uAtSecondGrating).^2)))});
        clim([0,ceil(max(max(abs(uAtSecondGrating).^2)))])
        title('Intensity at Second Grating Plane', 'FontWeight', 'normal')
        p1 = get(s1,'Position');
        p1(1) = 0.05; p1(3) = p1(3);
        set(s1,'Position',p1);

        % Show phase pattern at the second grating plane
        subplot(s2);
        imagesc(x0*1000,y0*1000,angle(uAtSecondGrating))
        colormap(s2,'jet'); axis image
        xlim(round(1000*[-wx/2,wx/2])); xticks(round(1000*[-wx/2,0,wx/2]))
        ylim(round(1000*[-wy/2,wy/2])); yticks(round(1000*[-wy/2,0,wy/2]))
        set(s2,'YDir','normal')
        colorbar(s2,'Ticks',[-pi,0,pi],'TickLabels',{'$-\pi$','0','$\pi$'}); clim([-pi,pi])
        title('Phase at Second Grating Plane', 'FontWeight', 'normal')
        p2 = get(s2,'Position');
        p2(1) = 0.05;
        set(s2,'Position',p2);
        
        % ------- [ Propagation to Lens | Free-Space + Lens Phase ] -------
        % Propagate from second grating to lens and apply thin lens phase.
        f = 20e-2; % focal length of the lens [m]
        zToLens = 5e-2; % distance from second grating to lens [m]
        lensPhase = exp(-i.*k.*r.^2./(2.*f));
        uAtLens = lensPhase .* ifft2(fft2(uAtSecondGrating).*exp(-i.*pi.*wl*zToLens.*(fx.^2+fy.^2)));
        
        % ------- [ Final Propagation | Detection Plane ] -------
        % Propagate from lens to detection plane (slightly after focal plane).
        zToDetectionPlane = f + 1e-2;
        uAtDetectionPlane = ifft2(fft2(uAtLens).*exp(-i.*pi.*wl*zToDetectionPlane.*(fx.^2+fy.^2)));
        
        % Show intensity pattern at the detection plane
        subplot(s3);
        imagesc(x0*1000,y0*1000,abs(uAtDetectionPlane).^2)
        colormap(s3,'hot'); axis image
        xlim([-0.4,0.4]); xticks([-0.4,0,0.4])
        ylim([-0.3,0.3]); yticks([-0.3,0,0.3])
        set(s3,'YDir','normal')
        clim([0, cbarMax]);  % FIXED limits
        colorbar(s3, 'Ticks', [0, cbarMax], 'TickLabels', {'0', num2str(cbarMax)});
        title('Intensity at Detection Plane ($z = f + 1$ cm)', 'FontWeight', 'normal')
        p3 = get(s3,'Position');
        p3(1) = 0.465; p3(3) = 0.45;
        set(s3,'Position',p3);
        
        % ----- Overlay intensity profile along horizontal midline -----
        hold(s3, 'on');
        
        % Extract intensity profile along y = 0
        intensityRow = abs(uAtDetectionPlane(round(Ny/2), :)).^2;
        
        % Normalize the profile to fit within the image height
        yOffset = -0.3;  % Position from bottom (adjust as needed)
        scaleFactor = 0.15;  % Height of the profile in mm (adjust as needed)
        profileScaled = yOffset + scaleFactor * (intensityRow / max(intensityRow));
        
        % Plot the profile
        plot(s3, x0*1000, profileScaled, 'w-', 'LineWidth', 2);
        l = line([-0.4 0.4],[0 0],'Color','w','LineWidth',1,'LineStyle','--');

        % Add vertical arrows on left and right sides
        hold(s3, 'on');
        
        % Left arrow:
        line(s3, [-0.35, -0.35], [0, -0.2], 'Color', 'w', 'LineWidth', 1, 'LineStyle','--');
        patch(s3, [-0.35-0.01, -0.35, -0.35+0.01], [-0.2, -0.2-0.02, -0.2], 'w');
        
        % Right arrow: 
        line(s3, [0.35, 0.35], [0, -0.2], 'Color', 'w', 'LineWidth', 1, 'LineStyle','--');
        patch(s3, [0.35-0.01, 0.35, 0.35+0.01], [-0.2, -0.2-0.02, -0.2], 'w');
        
        hold(s3, 'off');

        set(findall(gcf,'Type','text'), 'Interpreter','latex', 'FontSize',15)
        set(findall(gcf,'Type','axes'), 'FontSize',15, 'TickLabelInterpreter', 'latex')
        set(findall(gcf,'Type','colorbar'), 'FontSize',15, 'TickLabelInterpreter', 'latex')
        
        textTitle = ['$\theta = ' num2str(angularDeviationDeg(ii)) '^{\circ}$'];
        sgtitle(textTitle,'Interpreter','latex','FontSize',20)
        drawnow
        set(fig,'Units','normalized','Position',[0.07 0.11 0.88 0.77])
        
        % Capture frame for video
        frame = getframe(fig);
        writeVideo(videoWriter, frame);
        
        % Export selected frames at specific angles
        if ismember(angularDeviationDeg(ii), [0, 4.5, 9, 18])
            fileName = sprintf('detection_plane_m1=%d_theta=%.1fdeg.pdf', m1, angularDeviationDeg(ii));
            exportgraphics(fig, fullfile(dirImage, fileName));
        end
        % pause(1)
        cla(s1,s2,s3)
        
        % Extract on-axis intensity and maximum value
        onAxisIntensity(ii) = abs(uAtDetectionPlane(round(Ny/2),round(Nx/2))).^2;
        maxIntensity(ii) = max(max(abs(uAtDetectionPlane).^2));
        
    end

    close(videoWriter);
    close(fig);

    % Export maximum intensity vs angular deviation
    maxIntensityTable = table(angularDeviationDeg(:), maxIntensity(:), ...
        'VariableNames', {'AngularDeviation_deg', 'maxIntensity'});
    fileName = sprintf('max_intensity_m1=%d.txt', m1);
    writetable(maxIntensityTable, fullfile(dirData, fileName), 'Delimiter', '\t');
    
    % Export on-axis intensity vs angular deviation
    onAxisIntensityTable = table(angularDeviationDeg(:), onAxisIntensity(:), ...
        'VariableNames', {'AngularDeviation_deg', 'onAxisIntensity'});
    fileName = sprintf('on_axis_intensity_m1=%d.txt', m1);
    writetable(onAxisIntensityTable, fullfile(dirData, fileName), 'Delimiter', '\t');

    % Plot normalized on-axis intensity vs angular deviation
    figCurve = figure('Position', [161, 556, 560, 195]);
    plot(onAxisIntensityTable.AngularDeviation_deg, ...
         onAxisIntensityTable.onAxisIntensity / max(maxIntensity), '-o');
    xlabel('Angular Deviation (deg)');
    ylabel('On-axis Intensity');
    grid off;
    xlim([0, 18]);
    xticks(0:4.5:18);
    fileName = sprintf('on_axis_intensity_m1=%d.pdf', m1);
    exportgraphics(figCurve, fullfile(dirImage, fileName));
    close(figCurve);

end

%% ------- [ Combined Plot: Intensity vs Angle ] -------
% Plot normalized on-axis intensity vs angular deviation for all m1 values
figCombined1 = figure('Position', [117, 149, 1030, 418]);
hold on;
legendNames = {};

for m1 = [5, 10, 30, 50]
    % Read on-axis intensity table for current m1
    fileName = sprintf('on_axis_intensity_m1=%d.txt', m1);
    onAxisIntensityTable = readtable(fullfile(dirData, fileName), 'Delimiter', '\t');
    
    % Read maximum intensity table for normalization
    fileName = sprintf('max_intensity_m1=%d.txt', m1);
    maxIntensityTable = readtable(fullfile(dirData, fileName), 'Delimiter', '\t');

    % Normalize by global maximum intensity
    normalizedIntensity = onAxisIntensityTable.onAxisIntensity / cbarMax;
    
    % Plot with different styles for each m1 value
    if m1 == 5
        plot(onAxisIntensityTable.AngularDeviation_deg, normalizedIntensity, 'k-o');
    elseif m1 == 10
        plot(onAxisIntensityTable.AngularDeviation_deg, normalizedIntensity, 'r-pentagram', ...
             'MarkerFaceColor', 'r');
    elseif m1 == 30
        plot(onAxisIntensityTable.AngularDeviation_deg, normalizedIntensity, '-^', ...
             'Color', "#77AC30", 'MarkerFaceColor', "#77AC30");
    else  % m1 == 50
        plot(onAxisIntensityTable.AngularDeviation_deg, normalizedIntensity, 'b-v', ...
             'MarkerFaceColor', 'b');
    end
    legendNames{end+1} = sprintf('m_1 = %d', m1);
end

% Formatting and export
grid off;
xlim([0, 18]);
xticks(0:4.5:18);
legend(legendNames, 'Location', 'best');
xlabel('Angular Deviation (deg)');
ylabel('On-axis Intensity (normalized)');
box on;
exportgraphics(figCombined1, fullfile(dirImage, 'Combined_OnAxis_Intensity_vs_Angle.pdf'));

%% ------- [ Combined Plot: Angle vs Intensity (Swapped Axes) ] -------
% Plot angular deviation vs normalized on-axis intensity (axes swapped)
% Useful for reading angular tolerance at a given intensity level
close all;
figCombined2 = figure('Position', [618, 159, 314, 520]);
hold on;
legendNames = {};

% Custom color palette for different m1 values
RGB = [colororder('gem12'); colororder('dye')];
colorIndices = [4, 8, 11, 5];

for m1 = [5, 10, 30, 50]
    % Read data for current m1
    fileName = sprintf('on_axis_intensity_m1=%d.txt', m1);
    onAxisIntensityTable = readtable(fullfile(dirData, fileName), 'Delimiter', '\t');
    
    fileName = sprintf('max_intensity_m1=%d.txt', m1);
    maxIntensityTable = readtable(fullfile(dirData, fileName), 'Delimiter', '\t');

    % Swap axes: x = normalized intensity, y = angular deviation
    x = onAxisIntensityTable.onAxisIntensity / cbarMax;
    y = onAxisIntensityTable.AngularDeviation_deg;
    
    % Special angles to highlight with markers
    specialAngles = [0, 4.5, 9, 18];

    if m1 == 5
        color = RGB(colorIndices(1),:);
        plot(x, y,':', 'Color', color,'LineWidth',0.4)
        legendNames{end+1} = sprintf('m_1 = %d', m1); % Store legend entry
        for jj = 1:length(specialAngles)
        idx = find(y == specialAngles(jj));
            if ~isempty(idx)
                h = plot(x(idx), y(idx), 'o', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', color, ...
                    'MarkerSize', 7);
                set(get(get(h,'Annotation'),'LegendInformation'), ...
                'IconDisplayStyle','off')
            end
        end

    elseif m1 == 10
        color = RGB(colorIndices(2),:);
        plot(x, y, '-','Color', color,'LineWidth',1)
        legendNames{end+1} = sprintf('m_1 = %d', m1); % Store legend entry
        for jj = 1:length(specialAngles)
        idx = find(y == specialAngles(jj));
            if ~isempty(idx)
                h = plot(x(idx), y(idx), 's', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', color, ...
                    'MarkerSize', 7);
                set(get(get(h,'Annotation'),'LegendInformation'), ...
                'IconDisplayStyle','off')
            end
        end

    elseif m1 == 30
        color = RGB(colorIndices(3),:);
        plot(x, y, '.-','Color', color,'LineWidth',0.8)
        legendNames{end+1} = sprintf('m_1 = %d', m1); % Store legend entry
        for jj = 1:length(specialAngles)
        idx = find(y == specialAngles(jj));
            if ~isempty(idx)
                h = plot(x(idx), y(idx), '^', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', color, ...
                    'MarkerSize', 7);
                set(get(get(h,'Annotation'),'LegendInformation'), ...
                'IconDisplayStyle','off')
            end
        end

    else
        color = RGB(colorIndices(4),:);
        plot(x, y, '-.','Color', color,'LineWidth',0.6)
        legendNames{end+1} = sprintf('m_1 = %d', m1); % Store legend entry
        for jj = 1:length(specialAngles)
        idx = find(y == specialAngles(jj));
            if ~isempty(idx)
                h = plot(x(idx), y(idx), 'd', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', color, ...
                    'MarkerSize', 7);
                set(get(get(h,'Annotation'),'LegendInformation'), ...
                'IconDisplayStyle','off')
            end
        end
    end

    
end

% Formatting and export
grid off;
ylim([0, rad2deg(pi/m2)]);
yticks(0:4.5:18);
xlim([-0.05, 0.4]);
xticks(0:0.05:0.35);
grid on
legend(legendNames, 'Location', 'eastoutside');
xlabel('On-axis Intensity (normalized)');
ylabel('Angular Deviation (deg)');
box on;
exportgraphics(figCombined2, fullfile(dirImage, 'Combined_Angle_vs_OnAxis_Intensity.pdf'));
