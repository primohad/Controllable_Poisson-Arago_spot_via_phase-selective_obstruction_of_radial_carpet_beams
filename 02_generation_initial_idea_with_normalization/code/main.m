% ----------------------------------------------------------
% Generate longitudinal evolution data of field after amplitude grating
% with/without thin lens. Stores XZ/YZ maps, transverse planes, and movie.
%
% Output: ../data/ (XZ/YZ/XY matrices), ../image/ (PDF/PNG figures),
%         ../movie/ (AVI file)
%
% Notes: User selects simulation mode via menu.
%        Frames captured per z-slice for movie assembly.
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
% Ensures output directory exists; extend here if more outputs are added (e.g., data/, logs/).

dirImage = fullfile('..', 'image');
dirData  = fullfile('..', 'data');
dirMovie = fullfile('..', 'movie');

% Create directories safely
if ~exist(dirImage, 'dir')
    mkdir(dirImage);
end
if ~exist(dirData, 'dir')
    mkdir(dirData);
end
if ~exist(dirMovie, 'dir')
    mkdir(dirMovie);
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

dirPrefix = fullfile(dirImage,prefixName);
if ~exist(dirPrefix, 'dir')
    mkdir(dirPrefix);
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
frames = struct('cdata',[],'colormap',[]);
fig = figure(Position=[602.00 266.00 560.00 420.00]);

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
    
    onAxisIntensity = zeros(size(zPropagation)); % I(0,0; z)
    uXZPlane = zeros(length(x0),length(zPropagation)); % XZ slice (y=0) recorded along x for each z
    uYZPlane = zeros(length(y0), length(zPropagation));  % YZ slice (x=0) recorded along y for each z

    movieName = [sprintf('m1 = %d - m2 = %d', m1, m2) ...
            ' | ' prefixName ...
            ' | ' sprintf('f = %d.avi', f*100)];
    
    %% ------- [ Propagation Loop With Lens ] -------
    zTargets = [0 f-5e-2, f, f+5e-2 , f-10e-2];  % target zPropagation values
    
    % % First pass: find global maximum intensity and corresponding z position
    % fprintf('Finding global maximum intensity...\n');
    % globalMaxIntensity = 0;
    % zAtGlobalMax = 0;
    
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
    
    % Second pass: propagate with visualization and capture frames
    for ii = 1:length(zPropagation)
            
        uPropagated = ifft2(fft2(uAtLensPlane) .* exp(-i .* pi .* wl * zPropagation(ii) .* (fx.^2 + fy.^2)));
        
        % At z = 0 (first sample), enforce identity propagation (no numerical phase applied)
        if ii == 1
            uPropagated = uAtLensPlane;
        end
    
        % uPropagated = uPropagated ./ sqrt(globalMaxIntensity);
        uPropagated = uPropagated ./ sqrt(max(max(abs(uPropagated).^2)));
        
        % Intensity at current plane and on-axis sample (center pixel)
        I = abs(uPropagated).^2;
        onAxisIntensity(ii) = I(round(Ny/2), round(Nx/2));
        
        % Visualization (xy) for the current z-plane
        imagesc(x0*1000, y0*1000, I);
        colormap hot; axis image; set(gca, 'YDir', 'Normal');
        limit = 2.5; clim([0 1])
        ylim([-limit limit]); xlim([-limit limit]);
        
        % % Save a PNG snapshot to ../image
        % xticks([]); yticks([]);  % Hide ticks for clean frame capture
        % fileName = sprintf('frame%d.png', ii);
        % set(fig, "Position", [602.00 266.00 560.00 420.00]);
        % exportgraphics(fig, fullfile(dirImage, prefixName, fileName), Resolution=300);
        
        % Cosmetic axes updates (after saving the clean snapshot)
        colorbar;
        xticks(-limit:limit:limit);
        set(gca, 'YTick', get(gca, 'XTick'));
        title(['z = ' num2str(1000 * zPropagation(ii)) ' mm'], 'FontSize', 15);
        xlabel('x (mm)'); ylabel('y (mm)');
        
        % Capture frame for movie assembly
        frames(end+1) = getframe(fig);
        
        % Extract midline (y=0) profile for XZ map at this z
        uXZPlane(:, ii) = abs(uPropagated(round(Ny/2), :)).^2;
        % Extract midline (x=0) profile for YZ map at this z
        uYZPlane(:, ii) = abs(uPropagated(:, round(Nx/2))).^2;
    
        pause(1);
        
        % Keep the transverse field at the focal plane for later plotting
        if zPropagation(ii) == f
            uTransversePlane = uPropagated;
        end
    
        if ismember(zPropagation(ii), zTargets)
            intensityName = sprintf('m1 = %d - m2 = %d | %s | f = %d cm | zPropagation = %.2f cm | Intensity XYplane.txt', ...
                                    m1, m2, prefixName, f*100, 100*zPropagation(ii));
            writematrix(I, fullfile(dirData, intensityName));
        end
    end

%% ------- [ Propagation Without Lens ] -------
else
    
    movieName = [sprintf('m1 = %d - m2 = %d', m1, m2) ...
            ' | ' prefixName ...
            sprintf('.avi')];
    
    % z-sweep for free-space propagation (0 → 2.5 m), including a special plane at 1.5 m.
    zPropagation = unique(sort([linspace(0,2.5,(2^9-1)) 1.5 2.0 1.0]));

    onAxisIntensity = zeros(size(zPropagation)); % I(0,0; z) along optical axis
    uXZPlane = zeros(length(x0),length(zPropagation)); % XZ slice (y = 0) at each z
    uYZPlane = zeros(length(y0), length(zPropagation)); % YZ slice (x = 0) at each z


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
    
    % Second pass: propagate with visualization and capture frames
    for ii = 1:length(zPropagation)
        
        uPropagated = ifft2(fft2(uAtSecondGratingPlane) .* exp(-i .* pi .* wl * zPropagation(ii) .* (fx.^2 + fy.^2)));
        
        % At z = 0 (first sample), enforce identity propagation (no numerical phase)
        if ii == 1
            uPropagated = uAtSecondGratingPlane;
        end
    
        % uPropagated = uPropagated ./ sqrt(globalMaxIntensity);
        uPropagated = uPropagated ./ sqrt(max(max(abs(uPropagated).^2)));
        
        % Intensity at current plane and the on-axis sample (center pixel)
        I = abs(uPropagated).^2;
        onAxisIntensity(ii) = I(round(Ny/2), round(Nx/2));
        
        % Visualization (xy) for the current z-plane
        imagesc(x0*1000, y0*1000, I);
        clim([0 1]);
        colormap hot; axis image; set(gca, 'YDir', 'Normal');
        
        % % Save a PNG snapshot to ../image
        % xticks([]); yticks([]);
        % fileName = sprintf('frame%d.png', ii);
        % set(fig, "Position", [602.00 266.00 560.00 420.00]);
        % exportgraphics(fig, fullfile(dirImage, prefixName, fileName), Resolution=300);
        
        % Cosmetic axes updates (after saving the clean snapshot)
        colorbar;
        xticks((-winHalf:winHalf:winHalf) .* 1000);
        set(gca, 'YTick', get(gca, 'XTick'));
        title(['z = ' num2str(1000 * zPropagation(ii)) ' mm'], 'FontSize', 15);
        xlabel('x (mm)'); ylabel('y (mm)');
        
        % Capture and store frame for the movie
        frames(end+1) = getframe(fig);
        
        % Extract midline (y=0) profile for XZ map at this z
        uXZPlane(:, ii) = abs(uPropagated(round(Ny/2), :)).^2;
        % Extract midline (x=0) profile for YZ map at this z
        uYZPlane(:, ii) = abs(uPropagated(:, round(Nx/2))).^2;
    
        pause(1);
        
        % Keep the transverse field at z = 1.5 m for later plotting
        if zPropagation(ii) == 1.5
            uTransversePlane = uPropagated;
        end
    
        if ismember(zPropagation(ii), zTargets)
            intensityName = sprintf('m1 = %d - m2 = %d | %s | zPropagation = %.2f cm | Intensity XYplane.txt', ...
                                    m1, m2, prefixName, 100*zPropagation(ii));
            writematrix(I, fullfile(dirData, intensityName));
        end
    end

end

%% ------- [ XZ Visualization ] -------
%   - Top subplot: XZ map (longitudinal evolution) using the midline slice data uXZPlane.
%   - Bottom subplot: On-axis (x=0) intensity trace versus z.
% Notes:
%   - uXZPlane holds |U(x, y=0; z)|^2 along x for each z-sample (columns).
%   - Axis units: z in cm, x in mm.


figXZ = figure(Position=[729 259 581 407]);

s(1) = subplot(2,1,1);
imagesc (zPropagation*100, x0*1000, (uXZPlane)) % XZ intensity map: z [cm] vs x [mm]
if useLens % Zoom for "with lens" case (±2.5 mm)
    ylim([-2.5 2.5])
    yticks([-2 0 2])
else % Wider view for "without lens" case (±5.5 mm)
    ylim([-5.5 5.5])
    yticks([-5 0 5])
end

colormap hot; set(gca,"YDir","Normal");
colorbar;
xlabel('z (cm)'); ylabel('x (mm)')

% On-axis (x = 0) intensity trace vs z (uses the mid row of uXZPlane).
s(2) = subplot(2,1,2);
midRowIntensityXZ = (uXZPlane(round(Nx/2),:));
plot(zPropagation*100,midRowIntensityXZ)

% Align subplot widths (match the right edge/widths)
p1 = s(1).Position;
p2 = s(2).Position;

p1(3) = p2(3); % Make top subplot as wide as bottom subplot
s(1).Position = p1;

% Export XZ composite figure (saved to current working directory)
fileName = [sprintf('m1 = %d - m2 = %d', m1, m2) ...
            ' | ' prefixName ...
            ' | ' sprintf('XZ plane.pdf')];

exportgraphics (figXZ, fullfile(dirImage,fileName))

% Export raw XZ intensity data as a .txt matrix into ../data directory
%   - Each column corresponds to one z-sample (propagation distance).
%   - Each row corresponds to an x-coordinate (mm).
%   - This allows re-plotting or quantitative analysis later, independent of MATLAB figures.
XZintensityName = [sprintf('m1 = %d - m2 = %d', m1, m2) ...
            ' | ' prefixName ...
            ' | ' sprintf('XZ plane.txt')];
writematrix(uXZPlane,fullfile(dirData,XZintensityName));

%% ------- [ YZ Visualization ] -------
%   - Top subplot: YZ map (longitudinal evolution) using the midline slice data uYZPlane.
%   - Bottom subplot: On-axis (y=0) intensity trace versus z.
% Notes:
%   - uYZPlane holds |U(x=0, y; z)|^2 along y for each z-sample (columns).
%   - Axis units: z in cm, y in mm.

figYZ = figure(Position=[729 259 581 407]);

s(1) = subplot(2,1,1);
imagesc (zPropagation*100, y0*1000, (uYZPlane)) % XZ intensity map: z [cm] vs y [mm]
if useLens % Zoom for "with lens" case (±2.5 mm)
    ylim([-2.5 2.5])
    yticks([-2 0 2])
else % Wider view for "without lens" case (±5.5 mm)
    ylim([-5.5 5.5])
    yticks([-5 0 5])
end

colormap hot; set(gca,"YDir","Normal");
colorbar;
xlabel('z (cm)'); ylabel('y (mm)')

% On-axis (y = 0) intensity trace vs z (uses the mid row of uYZPlane).
s(2) = subplot(2,1,2);
midRowIntensityYZ = (uYZPlane(round(Ny/2),:));
plot(zPropagation*100,midRowIntensityYZ)

% Align subplot widths (match the right edge/widths)
p1 = s(1).Position;
p2 = s(2).Position;

p1(3) = p2(3); % Make top subplot as wide as bottom subplot
s(1).Position = p1;

% Export XZ composite figure (saved to current working directory)
fileName = [sprintf('m1 = %d - m2 = %d', m1, m2) ...
            ' | ' prefixName ...
            ' | ' sprintf('YZ plane.pdf')];

exportgraphics (figYZ, fullfile(dirImage,fileName))

% Export raw YZ intensity data as a .txt matrix into ../data directory
%   - Each column corresponds to one z-sample (propagation distance).
%   - Each row corresponds to a y-coordinate (mm).
%   - This allows re-plotting or quantitative analysis later, independent of MATLAB figures.
YZintensityName = [sprintf('m1 = %d - m2 = %d', m1, m2) ...
            ' | ' prefixName ...
            ' | ' sprintf('YZ plane.txt')];
writematrix(uYZPlane, fullfile(dirData,YZintensityName));

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
for ff = 2:length(frames)
    % convert the image to a frame
    frame = frames(ff) ;    
    writeVideo(writerObj, frame);
end

% Close the writer object (finalize the file)
close(writerObj);

%% ------- [ Transverse Plane Visualization ] -------
% Display transverse (xy) intensity distribution at selected z-plane
% Save figure to ../image with proper naming based on lens/no-lens branch
% uTransversePlane is propagated field stored at z = f (with lens) or z = 150 cm (without lens)

xyPlaneIntensity = abs(uTransversePlane).^2;

figTransverse = figure;

imagesc(x0*1000,y0*1000,xyPlaneIntensity);
axis image;
colormap hot;
colorbar
clim([0 round(max(max(xyPlaneIntensity)))])
set(gca,'YDir','Normal');
hold on

if useLens
    % Zoom into the focal spot region (±0.35 mm window)
    xlim([-0.2 0.2])
    ylim([-0.2 0.2])
    plot(x0*1000,(0.1)*xyPlaneIntensity(round(Ny/2),:)-0.1995,'w','LineWidth',2)
    fileName = [sprintf('m1 = %d - m2 = %d', m1, m2) ...
            ' | ' prefixName ...
            ' | ' sprintf('xy Plane at z = f.pdf')];
else
    % Coarser ticks for the no-lens propagation case
    xticks([-4 0 4])
    xticklabels([-4 0 4])
    yticks([-4 0 4])
    yticklabels([-4 0 4])
    xlim([-4.5 4.5]); ylim([-4.5 4.5])
    plot(x0*1000,2*xyPlaneIntensity(round(Ny/2),:)-4.494,'w','LineWidth',2)
    fileName = [sprintf('m1 = %d - m2 = %d', m1, m2) ...
            ' | ' prefixName ...
            ' | ' sprintf('xy Plane at z = 150 cm.pdf')];
end

xlabel('x (mm)'); ylabel('y (mm)')

% Export figure (saved under ../image)
exportgraphics(figTransverse,fullfile(dirImage,fileName));

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