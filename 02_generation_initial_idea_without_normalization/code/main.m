% ----------------------------------------------------------
% MATLAB R2025b
%
% Simulate z-evolution of radial carpet beam after amplitude grating.
%
% Workflow: Generate RCB -> MIS mask -> Amplitude grating ->
%           Propagate with/without lens -> Capture frames -> Assemble movie
%
% Output:
%   - ../movie/ (AVI movie)
%   - ../image/on_axis_intensity_curves_<mode>.pdf (dual-axis plot)
%
% Notes: User selects mode via menu.
%        Top subplot: beam image (3/4 height).
%        Middle: raw on-axis intensity.
%        Bottom: normalized on-axis intensity.
% ----------------------------------------------------------

%% ------- [ Environment Setup ] -------
% Clear workspace, close figures, compact display.

clc
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
% Define and create output directories

dirMovie = fullfile('..', 'movie');
dirImage = fullfile('..', 'image');

% Create directories safely
if ~exist(dirMovie, 'dir')
    mkdir(dirMovie);
end
if ~exist(dirImage, 'dir')
    mkdir(dirImage);
end

%% ------- [ Core Constants & Grid ] -------
% Define complex unit, grid sizes, and wavelength (SI units)

i = complex(0,1);
Nx = 2^11 + 1;  % Grid size along x
Ny = 2^11 + 1;  % Grid size along y
wl = 532e-9;    % Wavelength in meters (532 nm)

%% ------- [ User Selection: Lens Mode ] -------
% Display menu for simulation mode selection

choice = menu('Select simulation mode:', 'With Lens', 'Without Lens');
if choice == 0
    error('Selection cancelled by user.');
end
useLens = (choice == 1);     % Logical flag for downstream branching

%% ------- [ Simulation Parameters ] -------
% Define spatial window half-size based on selected mode

if useLens
    winHalf = 5e-3;          % [m] half-width for "with lens" scenario
    prefixName = 'with lens';
    disp('Simulation will include a lens in the propagation.');
else
    winHalf = 10e-3;         % [m] half-width for "without lens" scenario
    prefixName = 'without lens';
    disp('Simulation will proceed without a lens (free-space propagation).');
end

%% ------- [ Spatial Coordinates ] -------
% Build Cartesian grid and corresponding polar coordinates

x0 = linspace(-winHalf, winHalf, Nx);
y0 = linspace(-winHalf, winHalf, Ny);
[x, y] = meshgrid(x0, y0);
[tet, r] = cart2pol(x, y);

Lx = x0(end) - x0(1);        % [m] window width along x
Ly = y0(end) - y0(1);        % [m] window height along y
dx = Lx / Nx;                % [m] sample spacing x
dy = Ly / Ny;                % [m] sample spacing y

%% ------- [ Frequency Coordinates ] -------
% Frequency axes corresponding to FFT conventions

dfx = 1 / (Nx * dx);         % [1/m] fundamental frequency along x
dfy = 1 / (Ny * dy);         % [1/m] fundamental frequency along y
fx0 = (-Nx/2 : Nx/2 - 1) * dfx;   % centered frequency vector x
fy0 = (-Ny/2 : Ny/2 - 1) * dfy;   % centered frequency vector y
fx0 = fftshift(fx0);
fy0 = fftshift(fy0);
[fx, fy] = meshgrid(fx0, fy0);    % 2D frequency grids

%% ------- [ Field at Second Grating Plane ] -------
% Generate complex amplitude U(r,theta; z) at z = 2.3 m

k  = 2*pi/wl;                              % [1/m] wavenumber
zToSecondGrating = 2.3;                    % [m] propagation to the amplitude grating plane
Vp = pi/2;                                 % [rad] phase modulation amplitude (binary phase grating)
m1 = 10;                                   % [—] number of spokes in the phase grating (generator)

uToSecondGratingPlane = complexAmplitudeGenerator( ...
    m1, Vp, wl, zToSecondGrating, r, tet);  % Complex field at the second grating plane

% Normalize field magnitude to unity maximum (keeps dynamic range stable for thresholds).
uToSecondGratingPlane = uToSecondGratingPlane ./ max(max(sqrt(abs(uToSecondGratingPlane).^2)));

%% ------- [ MIS Annulus Mask ] -------
% Build annulus mask based on horizontal midline intensity profile
% Inner radius: point where intensity = 0.25 * first minimum
% Outer radius: first local minimum after the main peak

% Extract intensity along horizontal midline (y=0, right half)
intensity = abs(uToSecondGratingPlane).^2;
intensityHorizontalMidline = intensity(find(y0==0),find(x0==0):end);

% Smooth intensity profile using Gaussian filter to reduce noise
windowSize = 11;  % Odd number recommended for symmetry
smoothesIntensityHorizontalMidline = smoothdata(intensityHorizontalMidline, 'gaussian', windowSize);

% Detect peaks (maxima) in smoothed intensity profile
[pks, locsMax] = findpeaks(smoothesIntensityHorizontalMidline, x0(find(x0==0):end));

% Detect troughs (minima) by inverting the intensity signal
[troughs, locsMin] = findpeaks(-smoothesIntensityHorizontalMidline, x0(find(x0==0):end));
troughs = -troughs;  % Restore actual minimum values

% Extract first minimum and first maximum positions
firstMinValue = troughs(1);
firstMaxValue = pks(1);

% Find point before first maximum where intensity = 0.25 * firstMinValue
smoothesIntensityHorizontalMidlineBeforeFirstMax = ...
    smoothesIntensityHorizontalMidline(1:find(smoothesIntensityHorizontalMidline==firstMaxValue));

[~, closestIndex] = min(abs(smoothesIntensityHorizontalMidlineBeforeFirstMax - 0.25 * firstMinValue));
locInnerRadius = x0(find(x0 == 0) + closestIndex - 1); % Inner boundary of MIS
locOuterRadius = locsMin(1); % Outer boundary of MIS (first minimum)

% Create binary annulus mask (1 for MIS region, 0 elsewhere)
criterion = double((r <= locOuterRadius) & (r >= locInnerRadius));
uAfterMaskMIS = uToSecondGratingPlane .* criterion;  % Apply MIS ring mask

%% ------- [ Radial Amplitude Grating ] -------
% Build binary radial amplitude grating (m2 spokes) and apply

m2 = m1; % Number of spokes in the amplitude grating
secondGrating = (1/2) * (1 + sign(cos(m2 * tet))); % Binary radial amplitude mask (sector selector)

uAtSecondGratingPlane = uAfterMaskMIS .* secondGrating; % Field after amplitude grating

%% ------- [ Propagation Setup ] -------
% Prepare frame container and branch on lens/no-lens propagation

% Pre-allocate container for captured movie frames (to be filled in the loop)
frames = {};

fig = figure('Position', [13, 1, 754, 814], ...
             'Visible', 'off', ...
             'Color', 'white', ...
             'Resize', 'off', ...          % Prevent resizing to maintain aspect ratio
             'Units', 'pixels');           % Use pixels for consistent sizing
axTop = subplot(5, 1, 1:3);  % Top axis: beam intensity image
axMid = subplot(5, 1, 4);    % Middle axis: raw on-axis intensity vs z
axBot = subplot(5, 1, 5);    % Bottom axis: normalized on-axis intensity vs z

%% ------- [ Propagation With Lens ] -------
if useLens
    
    % --- Free-space propagation from amplitude grating to lens plane (zToLens) ---
    zToLens = 5e-2; % [m] distance from second grating to the lens
    uToLensPlane = ifft2(fft2(uAtSecondGratingPlane).*exp(-i .* pi .* wl * zToLens .* (fx.^2 + fy.^2)));
    
    % --- Thin lens phase transformation ---
    f = 20.*1e-2;  % [m] focal length of the lens (20 cm)
    lensPhase = exp(-i .* k .* r.^2 ./ (2.*f)); % Ideal thin-lens quadratic phase
    uAtLensPlane = lensPhase .* uToLensPlane; % Field immediately after the lens
    
    % Build a z-sweep tightly sampling around z = f (±5 cm window + endpoints)
    zPropagation = [linspace(0,f-(5e-2),2^8) linspace(f-(5e-2),f+(5e-2),2^8) f f-10e-2];
    zPropagation = sort(unique(zPropagation)); % Ensure monotonic, unique sampling

    movieName = [sprintf('m1 = %d - m2 = %d', m1, m2) ...
            ' | ' prefixName ...
            ' | ' sprintf('f = %d.avi', f*100)];
    
    %% ------- [ Propagation Loop With Lens ] -------
    % Target z positions for exporting selected frames
    zTargets = [0, f-5e-2, f, f+5e-2, f-10e-2];
    
    % % --- First pass: find global maximum intensity (commented out, kept for reference) ---
    % % This pass would find the maximum intensity across all z planes
    % % for normalization purposes.
    % fprintf('Finding global maximum intensity...\n');
    % globalMaxIntensity = 0;
    % zAtGlobalMax = 0;
    % 
    % for ii = 1:length(zPropagation)
    %     uPropagated = ifft2(fft2(uAtLensPlane) .* exp(-i .* pi .* wl * zPropagation(ii) .* (fx.^2 + fy.^2)));
    % 
    %     % At z = 0 (first sample), enforce identity propagation
    %     if ii == 1
    %         uPropagated = uAtLensPlane;
    %     end
    % 
    %     % Compute intensity at current plane
    %     I = abs(uPropagated).^2;
    %     currentMax = max(I(:));
    % 
    %     % Update global maximum if current plane has higher intensity
    %     if currentMax > globalMaxIntensity
    %         globalMaxIntensity = currentMax;
    %         zAtGlobalMax = zPropagation(ii);
    %     end
    % end
    % 
    % fprintf('Global maximum intensity: %.4f at z = %.2f mm\n', globalMaxIntensity, 1000 * zAtGlobalMax);
    
    % Preallocate on-axis intensity array
    onAxisIntensity = zeros(length(zPropagation), 1);

    % --- Initialize middle plot (raw on-axis intensity curve) ---
    % This shows the actual on-axis intensity (not normalized)
    hLineMid = plot(axMid, NaN, NaN, 'r-', 'LineWidth', 2);
    grid(axMid, 'on');
    xlim(axMid, [0, max(zPropagation) * 1000]);  % z-range in mm
    xticklabels('')  % Hide x-ticks to avoid clutter (x-label shown on bottom axis)
    ylim(axMid, [0, Inf]);

    % --- Initialize bottom plot (normalized on-axis intensity curve) ---
    % Intensity is normalized by the peak intensity of the current plane
    hLineBot = plot(axBot, NaN, NaN, 'r-', 'LineWidth', 2);
    xlabel(axBot, '$z$ (mm)', 'Interpreter', 'latex', 'FontSize', 12);
    grid(axBot, 'on');
    xlim(axBot, [0, max(zPropagation) * 1000]);  % z-range in mm
    ylim(axBot, [0, Inf]);

    % Second pass: propagate with visualization and capture frames
    for ii = 1:length(zPropagation)
            
        uPropagated = ifft2(fft2(uAtLensPlane) .* exp(-i .* pi .* wl * zPropagation(ii) .* (fx.^2 + fy.^2)));
        
        % At z = 0 (first sample), enforce identity propagation (no numerical phase applied)
        if ii == 1
            uPropagated = uAtLensPlane;
        end
    
        % uPropagated = uPropagated ./ sqrt(globalMaxIntensity);
        % uPropagated = uPropagated ./ sqrt(max(max(abs(uPropagated).^2)));
        
        % Intensity at current plane and on-axis sample (center pixel)
        I = abs(uPropagated).^2;

        % Extract intensity profile along y = 0
        intensityRow = I(round(Ny/2), :);

        % Store on-axis intensity
        onAxisIntensity (ii) = I(round(Ny/2), round(Nx/2));

        % Normalize on-axis intensity by the peak intensity of the current horizontal profile
        onAxisIntensityToMax(ii) = onAxisIntensity(ii) / max(intensityRow);
    
        % ----- Update top axis (beam image) -----
        axes(axTop);
        imagesc(axTop, x0*1000, y0*1000, I);
        colormap(axTop, 'hot');
        axis(axTop, 'image');
        set(axTop, 'YDir', 'normal');
        limit = 3;
        clim(axTop, [0, Inf]);
        xlim(axTop, [-limit, limit]);
        ylim(axTop, [-limit, limit]);
        hold(axTop, 'on');
        
        % Plot the profile
        yOffset = -limit + 0.05;
        scaleFactor = 1;
        profileScaled = yOffset + scaleFactor * (intensityRow / max(intensityRow));
        plot(axTop, x0*1000, profileScaled, 'w-', 'LineWidth', 1);
        line(axTop, [-limit, limit], [0, 0], 'Color', 'w', 'LineWidth', 0.5, 'LineStyle', ':');
        
        % Add vertical arrows on left and right sides
        % Left arrow
        line(axTop, [-(limit-0.5), -(limit-0.5)], [0, -2], 'Color', 'w', 'LineWidth', 0.5, 'LineStyle', '--');
        patch(axTop, [-(limit-0.5)-0.05, -(limit-0.5), -(limit-0.5)+0.05], ...
              [-2, -2-0.1, -2], 'w', 'EdgeColor', 'w', 'LineWidth', 1.5);
        
        % Right arrow
        line(axTop, [limit-0.5, limit-0.5], [0, -2], 'Color', 'w', 'LineWidth', 0.5, 'LineStyle', '--');
        patch(axTop, [(limit-0.5)-0.05, (limit-0.5), (limit-0.5)+0.05], ...
              [-2, -2-0.1, -2], 'w', 'EdgeColor', 'w', 'LineWidth', 1.5);
        
        hold(axTop, 'off');
        
        % Cosmetic updates for top axis
        colorbar(axTop);
        xticks(axTop, -limit:limit:limit);
        set(axTop, 'YTick', get(axTop, 'XTick'));
        title(axTop, ['$z = ' num2str(1000 * zPropagation(ii)) '$ mm']);
        xlabel(axTop, '$x$ (mm)');
        ylabel(axTop, '$y$ (mm)');
        
        % ----- Update middle axis (on-axis intensity curve) -----
        axes(axMid);
        set(hLineMid, 'XData', zPropagation(1:ii)*1000, ...
                   'YData', (onAxisIntensity(1:ii)));
        % ylim(axBot, [0, globalMaxIntensity]);
        
        % ----- Update bottom axis (on-axis intensity curve) -----
        axes(axBot);
        set(hLineBot, 'XData', zPropagation(1:ii)*1000, ...
                   'YData', (onAxisIntensityToMax(1:ii)));
        % ylim(axBot, [0, globalMaxIntensity]);

        % Apply LaTeX styling
        applyLatexToFigure(fig);

        % Adjust subplot positions to achieve desired layout (top large, middle/bottom small)
        axTop.Position = [0.13 0.50 0.69 0.47];
        axMid.Position = [0.13 0.25 0.78 0.12];
        axBot.Position = [0.13 0.06 0.78 0.12];
        
        % Ensure fixed figure size (re-apply to be safe)
        set(fig, 'Position', [13, 1, 754, 814]);
        frame = getframe(fig);
        frame.cdata = frame2im(frame);
        
        % Store first frame size and resize others
        if ii == 1
            targetSize = size(frame.cdata);
            targetSize = targetSize(1:2);
        end
        if size(frame.cdata,1) ~= targetSize(1) || size(frame.cdata,2) ~= targetSize(2)
            frame.cdata = imresize(frame.cdata, [targetSize(1), targetSize(2)]);
        end
        
        frames{end+1} = frame;
        
        pause(1);
    end

%% ------- [ Propagation Without Lens ] -------
else
    
    movieName = [sprintf('m1 = %d - m2 = %d', m1, m2) ...
            ' | ' prefixName ...
            sprintf('.avi')];
    
    % z-sweep for free-space propagation (0 → 2.5 m), including a special plane at 1.5 m.
    zPropagation = unique(sort([linspace(0,2.5,(2^9-1)) 1.5 2.0 1.0]));
    
    %% ------- [ Propagation Loop Without Lens ] -------
    zTargets = [0 1.0, 1.5, 2.0 2.5];
    % % First pass: find global maximum intensity and corresponding z position
    % fprintf('Finding global maximum intensity...\n');
    % globalMaxIntensity = 0;
    % zAtGlobalMax = 0;
    
    % for ii = 1:length(zPropagation)
    %     uPropagated = ifft2(fft2(uAtSecondGratingPlane) .* exp(-i .* pi .* wl * zPropagation(ii) .* (fx.^2 + fy.^2)));
    % 
    %     % At z = 0 (first sample), enforce identity propagation
    %     if ii == 1
    %         uPropagated = uAtSecondGratingPlane;
    %     end
    % 
    %     % Compute intensity at current plane
    %     I = abs(uPropagated).^2;
    %     currentMax = max(I(:));
    % 
    %     % Update global maximum if current plane has higher intensity
    %     if currentMax > globalMaxIntensity
    %         globalMaxIntensity = currentMax;
    %         zAtGlobalMax = zPropagation(ii);
    %     end
    % end
    % 
    % fprintf('Global maximum intensity: %.4f at z = %.2f mm\n', globalMaxIntensity, 1000 * zAtGlobalMax);
    
    % Preallocate on-axis intensity array
    onAxisIntensity = zeros(length(zPropagation), 1);

    % --- Initialize middle plot (raw on-axis intensity curve) ---
    % This shows the actual on-axis intensity (not normalized)
    hLineMid = plot(axMid, NaN, NaN, 'r-', 'LineWidth', 2);
    grid(axMid, 'on');
    xlim(axMid, [0, max(zPropagation) * 100]);  % z-range in mm
    xticklabels('')  % Hide x-ticks to avoid clutter (x-label shown on bottom axis)
    ylim(axMid, [0, Inf]);

    % --- Initialize bottom plot (normalized on-axis intensity curve) ---
    % Intensity is normalized by the peak intensity of the current plane
    hLineBot = plot(axBot, NaN, NaN, 'r-', 'LineWidth', 2);
    xlabel(axBot, '$z$ (cm)', 'Interpreter', 'latex', 'FontSize', 12);
    grid(axBot, 'on');
    xlim(axBot, [0, max(zPropagation) * 100]);  % z-range in mm
    ylim(axBot, [0, Inf]);

    % Second pass: propagate with visualization and capture frames
    for ii = 1:length(zPropagation)
        
        uPropagated = ifft2(fft2(uAtSecondGratingPlane) .* exp(-i .* pi .* wl * zPropagation(ii) .* (fx.^2 + fy.^2)));
        
        % At z = 0 (first sample), enforce identity propagation (no numerical phase)
        if ii == 1
            uPropagated = uAtSecondGratingPlane;
        end
            
        % Intensity at current plane and the on-axis sample (center pixel)
        I = abs(uPropagated).^2;
        
        % Extract intensity profile along y = 0
        intensityRow = I(round(Ny/2), :);

        % Store on-axis intensity
        onAxisIntensity (ii) = I(round(Ny/2), round(Nx/2));

        % Normalize on-axis intensity by the peak intensity of the current horizontal profile
        onAxisIntensityToMax(ii) = onAxisIntensity(ii) / max(intensityRow);

        % ----- Update top axis (beam image) -----
        axes(axTop);
        imagesc(axTop, x0*1000, y0*1000, I);
        colormap(axTop, 'hot');
        axis(axTop, 'image');
        set(axTop, 'YDir', 'normal');
        xticks((-winHalf:winHalf:winHalf) .* 1000);
        limit = winHalf * 1000;
        clim(axTop, [0, Inf]);
        hold(axTop, 'on');

        % Visualization (xy) for the current z-plane
        imagesc(x0*1000, y0*1000, I);
        clim([0 Inf]);
        colormap hot; axis image; set(gca, 'YDir', 'Normal');
        
        % Plot the profile
        yOffset = -limit + 0.05;
        scaleFactor = 5;
        profileScaled = yOffset + scaleFactor * (intensityRow);
        plot(axTop, x0*1000, profileScaled, 'w-', 'LineWidth', 1);
        line(axTop, [-limit, limit], [0, 0], 'Color', 'w', 'LineWidth', 0.5, 'LineStyle', ':');
        
        % Add vertical arrows on left and right sides
        % Left arrow
        line(axTop, [-(limit-2), -(limit-2)], [0, -8], 'Color', 'w', 'LineWidth', 0.5, 'LineStyle', '--');
        patch(axTop, [-(limit-2)-0.2, -(limit-2), -(limit-2)+0.2], ...
            [-8, -8-0.3, -8], 'w', 'EdgeColor', 'w', 'LineWidth', 1.5);

        % Right arrow
        line(axTop, [limit-2, limit-2], [0, -8], 'Color', 'w', 'LineWidth', 0.5, 'LineStyle', '--');
        patch(axTop, ...
            [(limit-2)-0.2, (limit-2), (limit-2)+0.2], ...   % base at x = limit-2
            [-8, -8-0.3, -8], ...                           % tip at y = -6.3 (pointing down)
          'w', 'EdgeColor', 'w', 'LineWidth', 1.5);
        hold(axTop, 'off');
        
        % Cosmetic updates for top axis
        colorbar(axTop);
        xticks((-winHalf:winHalf:winHalf) .* 1000);
        set(axTop, 'YTick', get(axTop, 'XTick'));
        title(axTop, ['$z = ' num2str(100 * zPropagation(ii)) '$ cm']);
        xlabel(axTop, '$x$ (mm)');
        ylabel(axTop, '$y$ (mm)');
        
        % ----- Update middle axis (on-axis intensity curve) -----
        axes(axMid);
        set(hLineMid, 'XData', zPropagation(1:ii)*100, ...
                   'YData', (onAxisIntensity(1:ii)));
        % ylim(axBot, [0, globalMaxIntensity]);
        
        % ----- Update bottom axis (on-axis intensity curve) -----
        axes(axBot);
        set(hLineBot, 'XData', zPropagation(1:ii)*100, ...
                   'YData', (onAxisIntensityToMax(1:ii)));
        % ylim(axBot, [0, globalMaxIntensity]);

        % Apply LaTeX styling
        applyLatexToFigure(fig);
        axTop.Position = [0.13 0.50 0.69 0.47];
        axMid.Position = [0.13 0.25 0.78 0.12];
        axBot.Position = [0.13 0.06 0.78 0.12];

        set(fig, 'Position', [13, 1, 754, 814]);
        frame = getframe(fig);
        frame.cdata = frame2im(frame);
        
        % Store first frame size and resize others
        if ii == 1
            targetSize = size(frame.cdata);
            targetSize = targetSize(1:2);
        end
        if size(frame.cdata,1) ~= targetSize(1) || size(frame.cdata,2) ~= targetSize(2)
            frame.cdata = imresize(frame.cdata, [targetSize(1), targetSize(2)]);
        end
        
        frames{end+1} = frame;
    
        pause(1);

    end

end

% Separate figure: raw intensity (left y) and normalized intensity (right y)
figCurves = figure('Color', 'white', 'Position', [26 76 1027 253]);
yyaxis left;
plot(zPropagation * (double(useLens) .* 1000 + double(~useLens) * 100), onAxisIntensity, 'b-', 'LineWidth', 2);
ylabel('On-axis Intensity (raw)', 'Interpreter', 'latex', 'FontSize', 12);

% z unit: mm for lens, cm for no lens
xlabel(sprintf('$z$ %s', char(useLens*'(mm)' + ~useLens*'(cm)')), 'Interpreter', 'latex', 'FontSize', 12);
grid on;

yyaxis right;
plot(zPropagation * (double(useLens) .* 1000 + double(~useLens) * 100), onAxisIntensityToMax, 'r-', 'LineWidth', 2);
ylabel('Normalized Intensity ($I/I_{\max}$)', 'Interpreter', 'latex', 'FontSize', 12);
ylim([0, 1.05]);

% Legend for dual-axis plot
legend({'Raw Intensity', 'Normalized Intensity'}, 'Location', 'northoutside', ...
       'Interpreter', 'latex', 'FontSize', 10);

applyLatexToFigure(figCurves, 12);
fileName = sprintf('on_axis_intensity_curves_%s.pdf', strrep(prefixName, ' ', '_'));
exportgraphics(figCurves, fullfile(dirImage, fileName));
% close(figCurves);

%% ------- [ Movie Assembly ] -------
% Collect frames stored in `movieFrames` and write to AVI file
% Uses movieName string set earlier

% Create the video writer object (AVI format)
writerObj = VideoWriter(fullfile(dirMovie,movieName));

% Set frame rate (10 fps)
writerObj.FrameRate = 10; % set the seconds per image

% Open the video writer
open(writerObj);

% Write each captured frame into the video file
for ff = 1:length(frames)
    % convert the image to a frame
    frame = frames{ff};
    writeVideo(writerObj, frame);
end

% Close the writer object (finalize the file)
close(writerObj);

%% ------- [ Local Functions ] -------

function complexAmplitude = complexAmplitudeGenerator(spokesNumber, PhaseModulationAmplitude, waveLength, propagationDistance, r, tet)
    i = complex(0, 1);
    k = 2 * pi / waveLength;
    rho = (k * r.^2) / (4 * propagationDistance);
    m = spokesNumber;
    Vp = PhaseModulationAmplitude;

    u = 0;
    for q = 1:2:20
        u = u + ((2/q) * sqrt(pi/2) * sin(Vp/2) * (-i)^((m/2-1)*q+1)) .* ...
            (besselj((q*m+1)/2,rho)+i*besselj((q*m-1)/2,rho)) .* ...
            cos(q*m.*tet);
    end
    u = (cos(Vp) + sqrt(rho) .* exp(i*rho) .* u);
    complexAmplitude = u;
end

function applyLatexToFigure(figHandle, fontSize)
    if nargin < 1, figHandle = gcf; end
    if nargin < 2, fontSize = 12; end

    % --- Axes properties ---
    allAxes = findall(figHandle, 'Type', 'axes');
    for ax = allAxes'
        set(ax, 'FontName', 'Times New Roman', 'FontSize', fontSize, ...
                'TickLabelInterpreter', 'latex', 'DefaultTextInterpreter', 'latex');
        
        % Axis labels
        xlabel(ax, get(get(ax,'XLabel'),'String'), 'Interpreter','latex','FontName','Times New Roman','FontSize',fontSize);
        ylabel(ax, get(get(ax,'YLabel'),'String'), 'Interpreter','latex','FontName','Times New Roman','FontSize',fontSize);
        
        % Legend
        leg = get(ax, 'Legend');
        if ~isempty(leg) && isvalid(leg)
            set(leg, 'Interpreter','latex','FontName','Times New Roman','FontSize',fontSize-2);
        end
        
        % Colorbar
        cb = get(ax, 'Colorbar');
        if ~isempty(cb) && isvalid(cb)
            set(cb, 'FontName','Times New Roman','FontSize',fontSize-2, 'TickLabelInterpreter','latex');
        end
    end

    % --- Titles: find all text objects that are in axes and near the top ---
    allText = findall(figHandle, 'Type', 'text');
    for txt = allText'
        parent = get(txt, 'Parent');
        % Only if parent is an axes
        if isa(parent, 'matlab.graphics.axis.Axes')
            pos = get(txt, 'Position');
            yLim = get(parent, 'YLim');
            % If the text is in the upper 20% of the axes, treat as title
            if pos(2) > yLim(2) - 0.2 * diff(yLim)
                set(txt, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize + 5);
            else
                set(txt, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fontSize);
            end
        end
    end

    % --- Global titles (suptitle, sgtitle) ---
    allTextGlobal = findall(figHandle, 'Type', 'text');
    for txt = allTextGlobal'
        tag = get(txt, 'Tag');
        if any(strcmp(tag, {'suptitle', 'sgtitle'}))
            set(txt, 'Interpreter','latex','FontName','Times New Roman','FontSize',fontSize+10);
        end
    end

    set(figHandle, 'Name', get(figHandle, 'Name'), 'NumberTitle', 'off');
end