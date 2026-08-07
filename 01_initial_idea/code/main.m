% ----------------------------------------------------------
% Generate multi-panel figure for radial carpet beam diffraction
% using annular MIS mask and binary radial amplitude grating.
%
% Output: Multi-panel PDF saved in ../image/
%
% Notes: r1/r2 estimated from mid-row intensity profile.
%        Assumes first minimum after thresholding exists.
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

% Create directories safely
if ~exist(dirImage, 'dir')
    mkdir(dirImage);
end

%% ------- [ Core Constants & Grid ] -------
% Define complex unit, grid sizes, and wavelength (SI units)

i = complex(0,1);
Nx = 2^11 + 1;  % Grid size along x
Ny = 2^11 + 1;  % Grid size along y
wl = 532e-9;    % Wavelength in meters (532 nm)

%% ------- [ Spatial Coordinates ] -------
% Build Cartesian grid and corresponding polar coordinates

winHalf = 5e-3; % Half-window size (meters). Axes later displayed in millimeters (×1000 factor).
x0 = linspace(-winHalf, winHalf, Nx);
y0 = linspace(-winHalf, winHalf, Ny);
[x, y] = meshgrid(x0, y0);
[tet, r] = cart2pol(x, y);

%% ------- [ Phase Grating ] -------
% Generate binary phase grating with m1 radial spokes
% Each spoke alternates phase between +gamma and -gamma
m1 = 10;          % Number of spokes (even number required for radial carpet beam)
gamma = pi/2;     % Phase modulation amplitude (radians)

% Create grating: sign(cos(m1*tet)) gives ±1, then multiply by gamma
t1 = exp(i * gamma * sign(cos(m1 * tet)));

% Visualize phase distribution of the grating
figure;
imagesc(x0*1000, y0*1000, angle(t1));
colormap jet;
axis image;
set(gca, 'YDir', 'normal');
clim([-pi/2 pi/2]);

% Configure colorbar with LaTeX tick labels
cb = colorbar;
cb.Ticks = [-pi/2 0 pi/2];
cb.TickLabels = {'-\frac{\pi}{2}', '0', '\frac{\pi}{2}'};
cb.TickLabelInterpreter = "latex";

% Set axes ticks in millimeters
set(gca, 'XTick', 1000 * (-winHalf:winHalf:winHalf));
set(gca, 'YTick', get(gca, 'XTick'));

fileName = fullfile(dirImage, sprintf('phaseRadialGrating_m1=%d_gamma=%.2fpi.pdf', m1, gamma/pi));
exportgraphics(gca, fileName)
close(gcf)

%% ------- [ Figure Layout ] -------
% The exported PDF is later split into 4 subfigures (each subfigure contains a pair of tiles).

fig = figure('color','white','Units','normalized','NumberTitle','off','Position',[0.1271 -0.0378 0.6951 0.4689]);
tiledlayout(2,4) % Layout: 2 rows × 4 columns = 8 panels (intensity/phase/masks/overlays).
% Panels map (left→right, top→bottom):
% t1: |U|^2 at second-grating plane (z = 2.3 m), normalized; axes in mm.
% t2: Wrapped phase ∠U at second-grating plane; [-π, π]; axes in mm.
% t3: |U|^2 after MIS annulus mask (r1 ≤ r ≤ r2); dashed circles = r1/r2 (mm).
% t4: ∠U after MIS annulus mask; dashed circles = r1/r2 (mm).
% t5: Binary radial amplitude grating (sector selector) with m2 spokes.
% t6: |U|^2 after applying the amplitude grating; dotted lines ≈ sector boundaries.
% t7: Phase overlay on π-coded blocked sectors (visualizing complement of selector).
% t8: ∠U after amplitude grating (final transmitted field).

%% ------- [ Field at Second Grating Plane ] -------
% Generate complex amplitude U(r,theta; z) at z = 2.3 m and extract
% MIS (Main Intensity Spike) and Tail regions within a single lobe

k  = 2*pi/wl;                              % [1/m] wavenumber
zToSecondGrating = 2.3;                    % [m] propagation to the amplitude grating plane
Vp = pi/2;                                 % [rad] phase modulation amplitude (binary phase grating)
m1 = 10;                                   % [—] number of spokes in the phase grating (generator)

uToSecondGratingPlane = complexAmplitudeGenerator( ...
    m1, Vp, wl, zToSecondGrating, r, tet);  % Complex field at the second grating plane

% Normalize field magnitude to unity maximum (keeps dynamic range stable for thresholds).
uToSecondGratingPlane = uToSecondGratingPlane ./ max(max(sqrt(abs(uToSecondGratingPlane).^2)));

% Intensity at second grating plane
t1 = nexttile;
imagesc(x0*1000,y0*1000,abs(uToSecondGratingPlane).^2);
axis image; colormap(t1,'hot');
set(gca,'YDir','normal'); clim([0 1]); colorbar;

% Phase at second grating plane
t2 = nexttile;
imagesc(x0*1000,y0*1000,angle(uToSecondGratingPlane));
axis image; colormap(t2,'jet');
set(gca,'YDir','normal'); clim([-pi pi]);
colorbar('Ticks',[-pi,0,pi],'TickLabels',{'-\pi','0','\pi'});

% Extract intensity along horizontal midline (y=0, right half)
intensity = abs(uToSecondGratingPlane).^2;
intensityHorizontalMidline = intensity(find(y0==0),find(x0==0):end);
xMidline = x0(x0 >= 0);

% Smooth intensity profile using Gaussian filter to reduce noise
windowSize = 11;  % Odd number recommended for symmetry
smoothesIntensityHorizontalMidline = smoothdata(intensityHorizontalMidline, 'gaussian', windowSize);

% Detect peaks (maxima) in smoothed intensity profile
[pks, locsMax] = findpeaks(smoothesIntensityHorizontalMidline, x0(find(x0==0):end));

% Detect troughs (minima) by inverting the intensity signal
[troughs, locsMin] = findpeaks(-smoothesIntensityHorizontalMidline, x0(find(x0==0):end));
troughs = -troughs;  % Restore actual minimum values

% First minimum intensity value (extracted from the original profile)
firstMinValue = intensityHorizontalMidline(smoothesIntensityHorizontalMidline == troughs(1));
firstMaxValue = pks(1);

% Find inner radius (where intensity equals firstMinValue)
idxFirstMax = find(smoothesIntensityHorizontalMidline == firstMaxValue, 1);
intensityBeforeMax = intensityHorizontalMidline(1:idxFirstMax);
[~, closestIdx] = min(abs(intensityBeforeMax - firstMinValue));
locInnerRadius = xMidline(closestIdx);
locOuterRadius = locsMin(1);  % First minimum after the peak

% Define single lobe angular region (between ±π/(2*m1))
maskLobe = (tet > -pi/(2*m1)) & (tet < pi/(2*m1));

% MIS (Main Intensity Spot) region: within one lobe, between inner/outer radii,
% and above intensity threshold
maskMIS = maskLobe ...
    .* ((x <= locOuterRadius) & (x >= locInnerRadius)) ...
    .* (intensity >= 0.15);
maskMIS = logical(maskMIS);

% Tail region: within same lobe, outside outer radius, and above intensity threshold
maskTail = maskLobe ...
    .* (x > locOuterRadius) ...
    .* (intensity >= 0.05);
maskTail = logical(maskTail);

% Visualize extracted MIS and Tail regions overlaid on intensity and phase
figure();

% Left: Intensity with region overlay
s1 = subplot(1,2,1);
imagesc(x0*1000,y0*1000,intensity .* (maskMIS|maskTail))
axis image; colormap(s1,'hot');
set(gca,'YDir','normal','XTick',[],'YTick',[])
clim([0 1])
hold on
contour(x0*1000, y0*1000, maskMIS|maskTail, [0.5 0.5], 'w-', 'LineWidth', 1);

% Right: Phase with region overlay
s2 = subplot(1,2,2);
imagesc(x0*1000,y0*1000,angle(uToSecondGratingPlane) .* (maskMIS|maskTail))
axis image; colormap(s2,'jet');
set(gca,'YDir','normal','XTick',[],'YTick',[])
clim([-pi pi])
hold on
contour(x0*1000, y0*1000, maskMIS|maskTail, [0.5 0.5], 'w-', 'LineWidth', 1);

% Export figure showing MIS and Tail regions
fileName = fullfile(dirImage, sprintf('MIS_Tail_regions_m1=%d.pdf', m1));
exportgraphics(gcf, fileName);
close(gcf);

%% ------- [ MIS Ring-Selection Criterion ] -------
% Build ring (annulus) mask based on middle-row intensity profile
% - midRowIntensity: 1D profile (right half) along the central row (y = 0).
% - r1: first x where intensity crosses the lower threshold (>= 0.1).
% - r2: first local minimum after the thresholded region (via findpeaks on negated signal).
% - windowThickness: radial ring thickness (r2 - r1).

midRowIntensity = abs(uToSecondGratingPlane(round(Ny/2),round(Nx/2):end)).^2;
midRowX = x0(round(Nx/2):end);

r1 = midRowX(find(midRowIntensity>=0.1,1));

% Use findpeaks on the negated signal to detect the first local minimum beyond threshold.
[~, r2] = findpeaks(-midRowIntensity.*double(midRowIntensity>=0.1),midRowX);
r2 = r2(1); % first detected minimum position (upper bound)

%{
Goal: define a circular annulus mask that keeps only the MIS band.
- Zero out pixels with radius > r2 (outside the ring).
- Zero out pixels with radius < r1 (inside the inner hole).
The passband is r1 <= r <= r2.
%}

criterion = double((r <= r2) & (r >= r1)); % [0/1] annulus mask (binary)
uAfterMaskMIS = uToSecondGratingPlane .* criterion; % Apply MIS ring mask

% Intensity after MIS annulus mask
t3 = nexttile;
imagesc(x0*1000,y0*1000,abs(uAfterMaskMIS).^2);
axis image; colormap(t3,'hot');
set(gca,'YDir','normal'); clim([0 1]); colorbar;

numPoints = 2^10;  % Number of points on the circular path

% Generate circular path coordinates centered at the origin
theta = linspace(0, 2*pi, numPoints);
x_circle_in = 1000 * r1 * cos(theta);
y_circle_in = 1000 * r1 * sin(theta);
x_circle_out = 1000 * r2 * cos(theta);
y_circle_out = 1000 * r2 * sin(theta);

% Overlay the circular path on the beam profile
hold on
plot(t3,x_circle_in, y_circle_in, '--w', 'LineWidth', 1);
plot(t3,x_circle_out, y_circle_out, '--w', 'LineWidth', 1);

% Phase after MIS annulus mask
t4 = nexttile;
imagesc(x0*1000,y0*1000,angle(uAfterMaskMIS));
axis image; colormap(t4,'jet');
set(gca,'YDir','normal'); clim([-pi pi]);
colorbar('Ticks',[-pi,0,pi],'TickLabels',{'-\pi','0','\pi'});

% Overlay the circular path on the beam profile
hold on
plot(t4,x_circle_in, y_circle_in, '--k', 'LineWidth', 1);
plot(t4,x_circle_out, y_circle_out, '--k', 'LineWidth', 1);

%% ------- [ Radial Amplitude Grating ] -------
% Brief: Build a binary radial amplitude grating (m2 spokes) and apply it to the
%        MIS-filtered field. Here m2 = m1 selects every other MIS with identical
%        phase structure (per the conceptual analysis).

m2 = m1; % Number of spokes in the amplitude grating
secondGrating = (1/2) * (1 + sign(cos(m2 * tet))); % Binary radial amplitude mask (sector selector)

% Binary amplitude grating (sector selector)
t5 = nexttile;
imagesc(x0*1000,y0*1000,abs(secondGrating).^2);
axis image; colormap(t5,'gray');
set(gca,'YDir','normal'); clim([0 1]); colorbar;

uAtSecondGratingPlane = uAfterMaskMIS .* secondGrating; % Field after amplitude grating

% Intensity after amplitude grating
t6 = nexttile;
imagesc(x0*1000,y0*1000,abs(uAtSecondGratingPlane).^2);
axis image; colormap(t6,'hot');
set(gca,'YDir','normal'); clim([0 1]); colorbar;

hold on
slope = zeros(1,m2);
slope(1) = pi/(2*m2);
y(1,:) = tan(slope(1)).*x0;
plot(t6,1000 * x0,1000 * y(1,:),':w','LineWidth',1)
for jj = 2:m2
    slope(jj) = slope(jj-1) + pi/m2;
    y(jj,:) = tan(slope(jj)) .* x0;
    plot(t6,1000 * x0, 1000 * y(jj,:),':w','LineWidth',1)
end
hold off
xlim([-winHalf winHalf]*1000); ylim([-winHalf winHalf]*1000)

% Phase overlay with blocked sectors
t7 = nexttile;
imagesc(x0*1000,y0*1000,pi .* double(~abs(secondGrating)))
axis image; colormap(t7,'jet');
set(gca,'YDir','normal'); clim([-pi pi]);
colorbar('Ticks',[-pi,0,pi],'TickLabels',{'-\pi','0','\pi'});
hold on
alpha = 0.7; % Set the transparency level
overlay = imagesc(x0*1000,y0*1000,angle(uAtSecondGratingPlane));
set(overlay, 'AlphaData', alpha);
hold off

% Phase after amplitude grating
t8 = nexttile;
imagesc(x0*1000,y0*1000,angle(uAtSecondGratingPlane))
axis image; colormap(t8,'jet');
set(gca,'YDir','normal'); clim([-pi pi]);
colorbar('Ticks',[-pi,0,pi],'TickLabels',{'-\pi','0','\pi'});

%% ------- [ Export Results ] -------
fileName = fullfile(dirImage,sprintf('RCB_MIS_transmission_m1=%d_m2=%d.pdf', m1, m2));
exportgraphics(fig, fileName)

%% ------- [ Local Functions ] -------

function complexAmplitude = complexAmplitudeGenerator(spokesNumber, PhaseModulationAmplitude, waveLength, propagationDistance, r, tet)
    % Generate complex field of radial carpet beam at propagation distance.

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