% ----------------------------------------------------------
% Identify in-phase/out-of-phase regions inside MIS (Main Intensity Spot)
% and Tail regions at the second grating plane for m1 = m2 = 10, 25, 50.
%
% Workflow:
%   1. Generate hybrid grating (m1) on SLM aperture
%   2. Propagate to second grating plane (z = 2.3 m)
%   3. Apply second binary radial amplitude grating (m2 = m1)
%   4. Define MIS and Tail regions based on intensity profile
%   5. Calculate phase statistics and power ratios
%   6. Export figures and numerical summaries%
% Output: ../image/ (MIS_tail figures, intensity profiles, phase summary)
%         ../data/  (Excel summaries per m1)
%
% Notes: m1 = m2 = [10, 25, 50].
%        MIS: inner region between r1 and r2 (first minimum).
%        Tail: region outside r2 within the same lobe.
%        In-phase = within mean phase ± STD.
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
% Define and create output directories

dirImage = fullfile('..', 'image');
dirData  = fullfile('..', 'data');

if ~exist(dirImage, 'dir'), mkdir(dirImage); end
if ~exist(dirData, 'dir'),  mkdir(dirData);  end

%% ------- [ Optical Constants ] -------
% Define complex unit, wavelength, and wavenumber

i = complex(0, 1);
wl = 532e-9;        % Wavelength in meters (532 nm)
k = 2*pi/wl;        % Wavenumber [1/m]

%% ------- [ Spatial Coordinates ] -------
% Build Cartesian grid and corresponding polar coordinates

Nx = 2^12 + 1;      % Grid size along x
Ny = 2^12 + 1;      % Grid size along y

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

wx = 1.6e-2;        % SLM width [m]
wy = 1.2e-2;        % SLM height [m]

uPlane = zeros(size(x));
maskSLM = (x > -wx/2 & x < wx/2) & (y > -wy/2 & y < wy/2);
uPlane(maskSLM) = 1;  % Unit amplitude over SLM aperture

%% ------- [ Grating Parameters ] -------
% Modulation parameters for hybrid grating

Vp = pi/2;          % Phase modulation amplitude [rad]
Va = 0.1;           % Amplitude visibility

%% ------- [ Main Loop: Sweep m1 = m2 ] -------
% Create complete figure for multi-panel display

figComplete = figure('Color', 'white', 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
tiledlayout(3, 2);

m1Values = [10, 25, 50];

for m1 = m1Values
    fprintf('Processing m1 = %d\n', m1);
    m2 = m1;  % m2 = m1 for this analysis

    % Create first hybrid radial grating (amplitude + phase)
    t1 = (0.5 * (1 + Va * sign(cos(m1 * tet)))) .* ...
         exp(1i * Vp * sign(cos(m1 * tet)));
    
    % Apply grating on SLM aperture
    uGratingOnSLM = t1 .* uPlane;
    
    % Propagate from SLM to second grating plane (z = 2.3 m)
    zPropagation = 2.3;  % [m]
    uPropagateBefore = ifft2(fft2(uGratingOnSLM) .* ...
        exp(-1i * pi * wl * zPropagation * (fx.^2 + fy.^2)));
    
    % Apply second binary radial amplitude grating (m2 = m1)
    t2 = 0.5 * (1 + sign(cos(m2 * tet)));
    uAfterSecondGrating = uPropagateBefore .* t2;
    
    % Normalize field magnitude
    uAfterSecondGrating = uAfterSecondGrating ./ ...
        sqrt(max(max(abs(uAfterSecondGrating).^2)));
    
    % Crop to square window matching SLM aperture
    [~, xIdxNeg] = min(abs(x0 - -wy/2));
    [~, xIdxPos] = min(abs(x0 - wy/2));
    [~, yIdxNeg] = min(abs(y0 - -wy/2));
    [~, yIdxPos] = min(abs(y0 - wy/2));
    
    uCropped = uAfterSecondGrating(yIdxNeg:yIdxPos, xIdxNeg:xIdxPos);
    x0Cropped = linspace(x0(xIdxNeg), x0(xIdxPos), size(uCropped, 2));
    y0Cropped = linspace(y0(yIdxNeg), y0(yIdxPos), size(uCropped, 1));
    
    % Compute intensity and phase of cropped field
    intensityCropped = abs(uCropped).^2;
    phaseCropped = angle(uCropped);
    
    % ----- Extract intensity profile along horizontal midline -----
    intensityMidline = intensityCropped(y0Cropped == 0, x0Cropped >= 0);
    xMidline = x0Cropped(x0Cropped >= 0);
    
    % Smooth intensity profile
    windowSize = round(Nx/100);  % Odd number for symmetry
    intensitySmooth = smoothdata(intensityMidline, 'gaussian', windowSize);
    
    % Find local maxima and minima
    [pks, locsMax] = findpeaks(intensitySmooth, xMidline, "MinPeakHeight", 0.1);
    [troughs, locsMin] = findpeaks(-intensitySmooth .* double(intensitySmooth>=0.1), xMidline);
    troughs = -troughs;  % Restore actual minimum values
    
    % Plot intensity profile with peaks and troughs
    figProfile = figure('Color', 'white');
    plot(xMidline * 1000, intensitySmooth, 'b', 'LineWidth', 1.5);
    hold on;
    plot(locsMax * 1000, pks, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    plot(locsMin * 1000, troughs, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
    plot(xMidline * 1000, intensityMidline, 'magenta', 'LineWidth', 0.5);
    xlabel('$x$ (mm)');
    ylabel('Intensity');
    title(sprintf('Intensity Profile - $m_1 =$ %d', m1));
    legend('Smoothed', 'Maxima', 'Minima', 'Original', 'Location', 'best');
    grid on;
    hold off;
    
    applyLatexToFigure(figProfile)
    
    % Export intensity profile figure
    fileName = sprintf('intensity_profile_m1_%d.pdf', m1);
    exportgraphics(figProfile, fullfile(dirImage, fileName));
    close(figProfile);
    
    % ----- Define MIS and Tail regions -----
    % First minimum intensity value (extracted from the original profile)
    firstMinValue = intensityMidline(intensitySmooth == troughs(1));
    firstMaxValue = pks(1);

    % Find inner radius (where intensity equals firstMinValue)
    idxFirstMax = find(intensitySmooth == firstMaxValue, 1);
    intensityBeforeMax = intensitySmooth(1:idxFirstMax);
    [~, closestIdx] = min(abs(intensityBeforeMax - firstMinValue));
    locInnerRadius = xMidline(closestIdx);
    locOuterRadius = locsMin(1);  % First minimum after the peak
    
    % Extract phase along midline between inner and outer radii
    idxInner = find(x0Cropped == locInnerRadius);
    idxOuter = find(x0Cropped == locOuterRadius);
    phaseMidline = phaseCropped(y0Cropped == 0, idxInner:idxOuter);
    
    % Define lobe mask (within ±π/m1)
    tetCropped = tet(yIdxNeg:yIdxPos, xIdxNeg:xIdxPos);
    maskLobe = (tetCropped >= -pi/(2*m1)) & (tetCropped <= pi/(2*m1));
    
    xCropped = x(yIdxNeg:yIdxPos, xIdxNeg:xIdxPos);
    
    % MIS region: lobe between inner/outer radii, intensity >= 0.15 (bright core)
    maskMIS = maskLobe & (xCropped <= locOuterRadius) & ...
              (xCropped >= locInnerRadius) & (intensityCropped >= 0.15);
    maskMIS = logical(maskMIS);
    
    % Tail region: lobe outside outer radius, intensity >= 0.05 (extended weak regions)
    maskTail = maskLobe & (xCropped > locOuterRadius) & (intensityCropped >= 0.05);
    maskTail = logical(maskTail);
    
    % ----- MIS statistics -----
    phaseMIS = phaseCropped(maskMIS);
    meanPhaseMIS = mean(phaseMIS) / pi;
    stdPhaseMIS = std(phaseMIS) / pi;
    
    phaseRangeMIS = pi * [meanPhaseMIS - stdPhaseMIS, meanPhaseMIS + stdPhaseMIS];
    inPhaseMIS = maskMIS & (phaseCropped > phaseRangeMIS(1)) & ...
                 (phaseCropped < phaseRangeMIS(2));
    outPhaseMIS = maskMIS & ~inPhaseMIS;
    
    inPhasePowerMIS = sum(intensityCropped(inPhaseMIS));
    outPhasePowerMIS = sum(intensityCropped(outPhaseMIS));
    percentageMIS = inPhasePowerMIS / outPhasePowerMIS;
    
    % ----- Tail statistics -----
    phaseTail = phaseCropped(maskTail);
    meanPhaseTail = mean(phaseTail) / pi;
    stdPhaseTail = std(phaseTail) / pi;
    
    phaseRangeTail = pi * [meanPhaseTail - stdPhaseTail, meanPhaseTail + stdPhaseTail];
    inPhaseTail = maskTail & (phaseCropped > phaseRangeTail(1)) & ...
                  (phaseCropped < phaseRangeTail(2));
    outPhaseTail = maskTail & ~inPhaseTail;
    
    inPhasePowerTail = sum(intensityCropped(inPhaseTail));
    outPhasePowerTail = sum(intensityCropped(outPhaseTail));
    percentageTail = inPhasePowerTail / outPhasePowerTail;
    
    % ----- Create MIS/Tail visualization figure -----
    figMISTail = figure('Color', 'white', 'Position', [356, 400, 1050, 234]);
    
    % Row 1: MIS region (columns 1-4)
    s(1) = subplot(2, 4, 1);
    imagesc(x0Cropped, y0Cropped, intensityCropped .* inPhaseMIS);
    axis equal; clim([0, 1]); colormap(s(1), 'hot'); xticks([]); yticks([]);
    applyLimOnXyMIS(maskMIS, x0Cropped, y0Cropped);
    title('Intensity (In-Phase MIS)', 'FontSize', 10);
    
    s(2) = subplot(2, 4, 5);
    imagesc(x0Cropped, y0Cropped, phaseCropped .* inPhaseMIS);
    axis equal; clim([-pi, pi]); colormap(s(2), 'jet'); xticks([]); yticks([]);
    applyLimOnXyMIS(maskMIS, x0Cropped, y0Cropped);
    title('Phase (In-Phase MIS)', 'FontSize', 10);

    s(3) = subplot(2, 4, 2);
    imagesc(x0Cropped, y0Cropped, intensityCropped .* outPhaseMIS);
    axis equal; clim([0, 1]); colormap(s(3), 'hot'); xticks([]); yticks([]);
    applyLimOnXyMIS(maskMIS, x0Cropped, y0Cropped);
    title('Intensity (Out-of-Phase MIS)', 'FontSize', 10);
    
    s(4) = subplot(2, 4, 6);
    imagesc(x0Cropped, y0Cropped, phaseCropped .* outPhaseMIS);
    axis equal; clim([-pi, pi]); colormap(s(4), 'jet'); xticks([]); yticks([]);
    applyLimOnXyMIS(maskMIS, x0Cropped, y0Cropped);
    title('Phase (Out-of-Phase MIS)', 'FontSize', 10);
    
    % Row 2: Tail region (columns 3-4 and 7-8)
    s(5) = subplot(2, 4, 3);
    imagesc(x0Cropped, y0Cropped, intensityCropped .* inPhaseTail);
    axis equal; clim([0, 1]); colormap(s(5), 'hot'); xticks([]); yticks([]);
    applyLimOnXyTail(maskTail, x0Cropped, y0Cropped);
    title('Intensity (In-Phase Tail)', 'FontSize', 10);
    
    s(6) = subplot(2, 4, 7);
    imagesc(x0Cropped, y0Cropped, phaseCropped .* inPhaseTail);
    axis equal; clim([-pi, pi]); colormap(s(6), 'jet'); xticks([]); yticks([]);
    applyLimOnXyTail(maskTail, x0Cropped, y0Cropped);
    title('Phase (In-Phase Tail)', 'FontSize', 10);
    
    s(7) = subplot(2, 4, 4);
    imagesc(x0Cropped, y0Cropped, intensityCropped .* outPhaseTail);
    axis equal; clim([0, 1]); colormap(s(7), 'hot'); xticks([]); yticks([]);
    applyLimOnXyTail(maskTail, x0Cropped, y0Cropped);
    title('Intensity (Out-of-Phase Tail)', 'FontSize', 10);
    
    s(8) = subplot(2, 4, 8);
    imagesc(x0Cropped, y0Cropped, phaseCropped .* outPhaseTail);
    axis equal; clim([-pi, pi]); colormap(s(8), 'jet'); xticks([]); yticks([]);
    applyLimOnXyTail(maskTail, x0Cropped, y0Cropped);
    title('Phase (Out-of-Phase Tail)', 'FontSize', 10);

    applyLatexToFigure(figMISTail)
    % Export MIS/Tail figure
    fileName = sprintf('MIS_tail_m1_%d.pdf', m1);
    exportgraphics(figMISTail, fullfile(dirImage, fileName));
    fileName = sprintf('MIS_tail_m1_%d.png', m1);
    exportgraphics(figMISTail, fullfile(dirImage, fileName), 'Resolution', 300);
    close(figMISTail);
    
    % ----- Add to complete figure -----
    figure(figComplete);
    
    % Intensity panel
    sIntensity = nexttile;
    imagesc(x0Cropped, y0Cropped, intensityCropped .* double(maskMIS + maskTail));
    axis image; clim([0, 1]); colormap(sIntensity, 'hot');
    applyLimOnXyLobe(maskLobe, x0Cropped, y0Cropped);
    hold on;
    contour(x0Cropped, y0Cropped, maskMIS, [0.5, 0.5], 'magenta-', 'LineWidth', 2);
    contour(x0Cropped, y0Cropped, maskTail, [0.5, 0.5], 'b-', 'LineWidth', 2);
    
    % Add gray overlay for in-phase regions
    rgbGray = zeros([size(maskMIS), 3]);
    rgbGray(:,:,1) = 0.9 * (inPhaseMIS | inPhaseTail);
    rgbGray(:,:,2) = 0.9 * (inPhaseMIS | inPhaseTail);
    rgbGray(:,:,3) = 0.9 * (inPhaseMIS | inPhaseTail);
    alphaCombined = (inPhaseMIS | inPhaseTail) * 0.50;
    h = imagesc(x0Cropped, y0Cropped, rgbGray);
    set(h, 'AlphaData', alphaCombined);
    title(sprintf('$m_1 =$ %d', m1), 'FontSize', 14);
    hold off;
    
    % Phase panel
    sPhase = nexttile;
    imagesc(x0Cropped, y0Cropped, phaseCropped .* double(maskMIS + maskTail));
    axis image; clim([-pi, pi]); colormap(sPhase, 'jet');
    applyLimOnXyLobe(maskLobe, x0Cropped, y0Cropped);
    hold on;
    contour(x0Cropped, y0Cropped, maskMIS, [0.5, 0.5], 'magenta-', 'LineWidth', 2);
    contour(x0Cropped, y0Cropped, maskTail, [0.5, 0.5], 'b-', 'LineWidth', 2);
    
    % Add dark overlay for in-phase regions
    rgbGray = zeros([size(maskMIS), 3]);
    rgbGray(:,:,1) = 0.1 * (inPhaseMIS | inPhaseTail);
    rgbGray(:,:,2) = 0.1 * (inPhaseMIS | inPhaseTail);
    rgbGray(:,:,3) = 0.1 * (inPhaseMIS | inPhaseTail);
    alphaCombined = (inPhaseMIS | inPhaseTail) * 0.60;
    h = imagesc(x0Cropped, y0Cropped, rgbGray);
    set(h, 'AlphaData', alphaCombined);
    hold off;

    % ----- Save numerical results -----
    resultsTable = table(["m1"; "Mean Phase MIS"; "Std Phase MIS"; ...
                         "Power In-Phase MIS"; "Power Out-of-Phase MIS"; "Percentage MIS"; ...
                         "Mean Phase Tail"; "Std Phase Tail"; ...
                         "Power In-Phase Tail"; "Power Out-of-Phase Tail"; "Percentage Tail"], ...
                        [m1; meanPhaseMIS; stdPhaseMIS; ...
                         inPhasePowerMIS; outPhasePowerMIS; percentageMIS; ...
                         meanPhaseTail; stdPhaseTail; ...
                         inPhasePowerTail; outPhasePowerTail; percentageTail]);
    
    excelFile = fullfile(dirData, sprintf('results_m1_%d.xlsx', m1));
    writetable(resultsTable, excelFile, 'Sheet', 'Results', 'WriteVariableNames', false);
    fprintf('  Saved results to: %s\n', excelFile);
end

applyLatexToFigure(figComplete)

% Export complete figure
fileName = 'MIS_tail_Complete.pdf';
exportgraphics(figComplete, fullfile(dirImage, fileName));
fprintf('Complete figure saved to: %s\n', fullfile(dirImage, fileName));

%% ------- [ Phase Summary Plot ] -------
% Read Excel summaries and plot mean phase ± STD for MIS and Tail

fileList = {'results_m1_10.xlsx', 'results_m1_25.xlsx', 'results_m1_50.xlsx'};

% Initialize arrays
m1Values = [];
meanPhaseMIS = []; stdPhaseMIS = [];
inPhasePowerMIS = []; outPhasePowerMIS = [];
meanPhaseTail = []; stdPhaseTail = [];
inPhasePowerTail = []; outPhasePowerTail = [];

% Extract data from each file
for i = 1:length(fileList)
    fileName = fullfile(dirData, fileList{i});
    data = readtable(fileName, 'ReadVariableNames', false);
    
    m1Values(end+1) = data{1, 2};
    meanPhaseMIS(end+1) = data{2, 2};
    stdPhaseMIS(end+1) = data{3, 2};
    inPhasePowerMIS(end+1) = data{4, 2};
    outPhasePowerMIS(end+1) = data{5, 2};
    meanPhaseTail(end+1) = data{7, 2};
    stdPhaseTail(end+1) = data{8, 2};
    inPhasePowerTail(end+1) = data{9, 2};
    outPhasePowerTail(end+1) = data{10, 2};
end

% Create summary plot
figSummary = figure('Color', 'white', 'Units', 'pixels', 'Position', [141, 360, 918, 388]);
hold on;

% Plot MIS (red circles)
errorbar(m1Values, meanPhaseMIS, stdPhaseMIS, 'ro', 'MarkerFaceColor', 'r', ...
         'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'MIS (Red)');

% Plot Tail (blue squares)
errorbar(m1Values, meanPhaseTail, stdPhaseTail, 'bs', 'MarkerFaceColor', 'b', ...
         'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Tail (Blue)');

% Add shaded regions
for ii = 1:length(m1Values)
    xRange = [m1Values(ii) - 1, m1Values(ii) + 1];
    
    if inPhasePowerMIS(ii) > inPhasePowerTail(ii)
        if meanPhaseMIS(ii) > meanPhaseTail(ii)
            yStart = meanPhaseMIS(ii) + stdPhaseMIS(ii);
            yEnd = meanPhaseMIS(ii) + stdPhaseMIS(ii) - 1;
        else
            yStart = meanPhaseMIS(ii) - stdPhaseMIS(ii);
            yEnd = meanPhaseMIS(ii) - stdPhaseMIS(ii) + 1;
        end
        fill([xRange(1), xRange(2), xRange(2), xRange(1)], ...
             [yStart, yStart, yEnd, yEnd], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    else
        if meanPhaseTail(ii) > meanPhaseMIS(ii)
            yStart = meanPhaseTail(ii) + stdPhaseTail(ii);
            yEnd = meanPhaseTail(ii) + stdPhaseTail(ii) - 1;
        else
            yStart = meanPhaseTail(ii) - stdPhaseTail(ii);
            yEnd = meanPhaseTail(ii) - stdPhaseTail(ii) + 1;
        end
        fill([xRange(1), xRange(2), xRange(2), xRange(1)], ...
             [yStart, yStart, yEnd, yEnd], ...
             'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end

% Labels and formatting
xlabel('$m_1$');
ylabel('Phase ($\pi$ units)');
xticks([10, 25, 50]);
legend('Location', 'eastoutside');
box on;
grid off;
hold off;

applyLatexToFigure(figSummary)

% Export summary figure
fileName = 'Mean_Phase_and_STD_MIS_Tail.pdf';
exportgraphics(figSummary, fullfile(dirImage, fileName));
fprintf('Summary figure saved to: %s\n', fullfile(dirImage, fileName));

%% ------- [ Generate LaTeX Table ] -------
% Read data for all m1 values (values are already in π units)

m1Values = [10, 25, 50];
MIS_meanPhase = zeros(1, 3);
MIS_stdPhase = zeros(1, 3);
MIS_pIn = zeros(1, 3);
MIS_pOut = zeros(1, 3);
Tail_meanPhase = zeros(1, 3);
Tail_stdPhase = zeros(1, 3);
Tail_pIn = zeros(1, 3);
Tail_pOut = zeros(1, 3);

for i = 1:length(m1Values)
    data = readtable(fullfile(dirData, sprintf('results_m1_%d.xlsx', m1Values(i))), ...
                     'ReadVariableNames', false);
    MIS_meanPhase(i) = data{2, 2};   % Already in π units
    MIS_stdPhase(i) = data{3, 2};    % Already in π units
    MIS_pIn(i) = data{4, 2};
    MIS_pOut(i) = data{5, 2};
    Tail_meanPhase(i) = data{7, 2};  % Already in π units
    Tail_stdPhase(i) = data{8, 2};   % Already in π units
    Tail_pIn(i) = data{9, 2};
    Tail_pOut(i) = data{10, 2};
end

% Calculate coherence percentages
MIS_percentage = 100 * MIS_pIn ./ (MIS_pIn + MIS_pOut);
Tail_percentage = 100 * Tail_pIn ./ (Tail_pIn + Tail_pOut);

filePath = fullfile(dirData, 'Phase_Power_Table.tex');
if exist(filePath, 'file')
    delete(filePath);
    fprintf('File deleted: %s\n', filePath);
end

% Open LaTeX file for writing
fid = fopen(fullfile(dirData, 'Phase_Power_Table.tex'), 'w');

% Write table header
fprintf(fid, '\\begin{table}[h]\n');
fprintf(fid, '\t\\caption{Phase and power distribution analysis for the in-phase and out-of-phase regions of the MIS and tail at the second grating plane for different values of $m_1$}\\label{tab1}\n');
fprintf(fid, '\t\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lcccccc}\n');
fprintf(fid, '\t\t\\toprule%%\n');
fprintf(fid, '\t\t& \\multicolumn{3}{@{}c@{}}{MIS} & \\multicolumn{3}{@{}c@{}}{Tail} \\\\\\cmidrule{2-4}\\cmidrule{5-7}%%\n');
fprintf(fid, '\t\t$m_1$ & $10$ & $25$ & $50$ & $10$ & $25$ & $50$ \\\\\n');
fprintf(fid, '\t\t\\midrule\n');

% Write Mean Phase row (values already in π units)
fprintf(fid, '\t\t$\\bar{\\phi}$  & $%.2f~\\pi$ & $%.2f~\\pi$ & $%.2f~\\pi$ & $%.2f~\\pi$ & $%.2f~\\pi$ & $%.2f~\\pi$\\\\\n', ...
        MIS_meanPhase(1), MIS_meanPhase(2), MIS_meanPhase(3), ...
        Tail_meanPhase(1), Tail_meanPhase(2), Tail_meanPhase(3));

% Write Std Phase row (values already in π units)
fprintf(fid, '\t\t$\\sigma_{\\phi}$  & $%.2f~\\pi$ & $%.2f~\\pi$ & $%.2f~\\pi$ & $%.2f~\\pi$ & $%.2f~\\pi$ & $%.2f~\\pi$\\\\\n', ...
        MIS_stdPhase(1), MIS_stdPhase(2), MIS_stdPhase(3), ...
        Tail_stdPhase(1), Tail_stdPhase(2), Tail_stdPhase(3));

% Write In-Phase Power row
fprintf(fid, '\t\t$P_{\\mathrm{in-phase}}$  & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$\\\\\n', ...
        MIS_pIn(1), MIS_pIn(2), MIS_pIn(3), ...
        Tail_pIn(1), Tail_pIn(2), Tail_pIn(3));

% Write Out-of-Phase Power row
fprintf(fid, '\t\t$P_{\\mathrm{out-of-phase}}$  & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$\\\\\n', ...
        MIS_pOut(1), MIS_pOut(2), MIS_pOut(3), ...
        Tail_pOut(1), Tail_pOut(2), Tail_pOut(3));

% Write Coherence Percentage row
fprintf(fid, '\t\t$100\\times\\frac{P_{\\mathrm{in-phase}}}{P_{\\mathrm{in-phase}}+P_{\\mathrm{out-of-phase}}}$  & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$\\\\\n', ...
        MIS_percentage(1), MIS_percentage(2), MIS_percentage(3), ...
        Tail_percentage(1), Tail_percentage(2), Tail_percentage(3));

% Write table footer
fprintf(fid, '\t\t\\bottomrule\n');
fprintf(fid, '\t\\end{tabular*}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid);
fprintf('LaTeX table saved to: %s\n', fullfile(dirData, 'Phase_Power_Table.tex'));

%% ------- [ Local Functions ] -------

function [] = applyLimOnXyMIS(maskMIS, x0Square, y0Square)
    % Restrict axes to bounding box of MIS mask
    [rowIdx, colIdx] = find(maskMIS);
    
    xMin = min(colIdx); xMax = max(colIdx);
    yMin = min(rowIdx); yMax = max(rowIdx);
    
    xlim(round([x0Square(xMin), x0Square(xMax)], 4));
    xticks(round([x0Square(xMin), x0Square(xMax)], 4));
    xticklabels(round([x0Square(xMin), x0Square(xMax)], 4) * 1000);
    ylim(round(y0Square(yMax), 4) * [-1 1]);
    yticks(round(y0Square(yMax), 4) * [-1 1]);
    yticklabels(round(y0Square(yMax), 4) * [-1 1] * 1000);
    set(gca, "YDir", "Normal");
end

function [] = applyLimOnXyTail(maskTail, x0Square, y0Square)
    % Restrict axes to bounding box of Tail mask
    [rowIdx, colIdx] = find(maskTail);
    
    xMin = min(colIdx); xMax = max(colIdx);
    yMin = min(rowIdx); yMax = max(rowIdx);
    
    xlim(round([x0Square(xMin), x0Square(xMax)], 4));
    xticks(round([x0Square(xMin), x0Square(xMax)], 4));
    xticklabels(round([x0Square(xMin), x0Square(xMax)], 4) * 1000);
    ylim(round(y0Square(yMax), 4) * [-1 1]);
    yticks(round(y0Square(yMax), 4) * [-1 1]);
    yticklabels(round(y0Square(yMax), 4) * [-1 1] * 1000);
    set(gca, "YDir", "Normal");
end

function [] = applyLimOnXyLobe(maskLobe, x0Square, y0Square)
    % Restrict axes to bounding box of Lobe mask
    [rowIdx, colIdx] = find(maskLobe);
    
    xMin = min(colIdx); xMax = max(colIdx);
    yMin = min(rowIdx); yMax = max(rowIdx);
    
    xlim(round([x0Square(xMin), x0Square(xMax)], 4));
    xticks(round([x0Square(xMin), x0Square(xMax)], 4));
    xticklabels(round([x0Square(xMin), x0Square(xMax)], 4) * 1000);
    ylim(round(y0Square(yMax), 4) * [-1 1]);
    yticks(round(y0Square(yMax), 4) * [-1 1]);
    yticklabels(round(y0Square(yMax), 4) * [-1 1] * 1000);
    set(gca, "YDir", "Normal");
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
    %
    % Notes:
    %   - Converts ALL text elements: axes labels, titles, legends, colorbars,
    %     text annotations, subplot titles, etc.
    %   - Sets font to Times New Roman
    %   - Sets interpreter to LaTeX for all text
    
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
            set(cb, 'FontName', 'Times New Roman', 'FontSize', fontSize - 2);
            set(cb, 'TickLabelInterpreter', 'latex');
            % Colorbar label
            cbLabel = get(cb, 'Label');
            if ~isempty(cbLabel) && isvalid(cbLabel)
                set(cbLabel, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize);
            end
        end
    end
    
    % ----- Text Annotations (if any) -----
    allText = findall(figHandle, 'Type', 'text');
    for i = 1:length(allText)
        set(allText(i), 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize);
    end
    
    % ----- Suptitle (if exists) -----
    allSuptitles = findall(figHandle, 'Type', 'text', 'Tag', 'suptitle');
    for i = 1:length(allSuptitles)
        set(allSuptitles(i), 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize + 4);
    end
    
    % ----- Figure Title (Name) -----
    set(figHandle, 'Name', get(figHandle, 'Name'), 'NumberTitle', 'off');
    
end