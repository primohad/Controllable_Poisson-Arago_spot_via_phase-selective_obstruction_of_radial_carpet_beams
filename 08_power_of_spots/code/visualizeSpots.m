% ----------------------------------------------------------
% Visualize precomputed central-spot and re-simulated full-field
% intensity patterns for m1 = [1,3,5]*m2.
%
% Workflow: Load spotU matrices -> Normalize to reference m1=m2 ->
%           Display central spots -> Re-simulate full fields ->
%           Normalize and display whole patterns
%
% Input: ../data/spotIntensity_m2=<m2>_m1=<m1>.txt
% Output: ../image/centralSpots_m2=<m2>.pdf/png
%         ../image/wholePattern_m2=<m2>.pdf
%
% Notes: m2 = 10.
%        Normalization relative to m1 = m2.
% ----------------------------------------------------------

% MATLAB R2025b

%% ------- [ Environment Setup ] -------
% Clear workspace, close figures, compact display

clc;
clear; close all; format compact; format bank;

%% ------- [ Output Paths ] -------
% Define and create output directories

dirImage  = fullfile('..', 'image');
dirData   = fullfile('..', 'data');

% Create directories safely
if ~exist(dirImage, 'dir'),  mkdir(dirImage);  end
if ~exist(dirData, 'dir'),   mkdir(dirData);   end

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

% Axis limits for zoomed views (mm)
xLim = 0.4/4; % [mm]
yLim = 0.3/3; % [mm]

%% ------- [ Visualization of Central Spot Patterns Only ] -------
% Loads and displays precomputed central spot patterns (spotIntensity) from file.
% Each image shows the normalized central region of intensity for m1 = [1×m2, 3×m2, 5×m2].

m2Values = 10;

for m2Idx = 1:length(m2Values)
    m2 = m2Values(m2Idx);
    fprintf('m2 = %d\n', m2);

    fig = figure('Position', [100, 100, 1000, 300]);

    % Load base matrix for normalization (m1 = 1×m2)
    baseM1 = m2 * 1;
    baseFilename = fullfile(dirData, sprintf('spotIntensity_m2=%d_m1=%d.txt', m2, baseM1));

    if isfile(baseFilename)
        baseSpotIntensity = readmatrix(baseFilename);
        normFactor = max(baseSpotIntensity(:));
    else
        warning('Base file not found: %s', baseFilename);
        normFactor = 1;
    end

    for m1Idx = 1:3
        ratio = 2 * m1Idx - 1;  % 1, 3, 5
        m1 = m2 * ratio;
        fprintf('\tm1 = %d * m2\n', ratio);

        filename = fullfile(dirData, sprintf('spotIntensity_m2=%d_m1=%d.txt', m2, m1));

        if isfile(filename)
            spotIntensity = readmatrix(filename);
            spotIntensity = spotIntensity / normFactor;
        else
            warning('File not found: %s', filename);
            continue;
        end

        subplot(1, 3, m1Idx);
        imagesc(x0*1e3, y0*1e3, spotIntensity);
        axis image;
        set(gca, 'YDir', 'normal', 'FontName', 'Times New Roman', 'FontSize', 10);
        colormap hot;
        xlim([-xLim, xLim]); ylim([-yLim, yLim]);
        title(sprintf('$m_1 = %d$', m1), 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 10);
        clim([0 1]);
        cb = colorbar;
        xticks(-xLim:xLim:xLim)
        yticks(-yLim:yLim:yLim)
    end

    sgtitle(sprintf('Normalized Central Spots for $m_2 = %d$', m2), ...
            'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 14);

    % Save and close figure
    exportgraphics(fig, fullfile(dirImage, sprintf('centralSpots_m2=%d.pdf', m2)));
    exportgraphics(fig, fullfile(dirImage, sprintf('centralSpots_m2=%d.png', m2)), 'Resolution', 300);
    close(fig);
end

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

%% ------- [ Full Field Simulation and Visualization ] -------
% Performs complete simulation from two radial gratings to image plane.
% Computes intensity profiles including central peak and surrounding rings.
% Each image shows the full simulated beam intensity pattern for m1 = [1×m2, 3×m2, 5×m2].

m2Values = 10;

for m2Idx = 1:length(m2Values)
    
    m2 = m2Values(m2Idx);
    fprintf('m2 = %d\n', m2);

    m1Values = m2 * [1 3 5];
    intensityCell = cell(1, length(m1Values));

    for m1Idx = 1:length(m1Values)
    
        m1 = m1Values(m1Idx);
        fprintf('\tm1 = %d\n', m1);
        
        % Create first hybrid radial grating (amplitude + phase)
        grating1 = (0.5 * (1 + Va * sign(cos(m1 * tet)))) .* ...
                   exp(1i * Vp * sign(cos(m1 * tet)));
        
        % Propagate from SLM to second grating plane (z = 2.3 m)
        zToSecondGrating = 2.3;  % [m]
        uAtSecondGrating = ifft2(fft2(grating1 .* uPlane) .* ...
            exp(-1i * pi * lambda * zToSecondGrating * (fx.^2 + fy.^2)));
     
        % Apply second binary radial amplitude grating (m2)
        grating2 = 0.5 * (1 + sign(cos(m2 * tet)));
        uAfterSecondGrating = uAtSecondGrating .* grating2;
        
        % Propagate from second grating to lens plane (z = 5 cm)
        zToLens = 5e-2;  % [m]
        uAtLens = ifft2(fft2(uAfterSecondGrating) .* ...
            exp(-1i * pi * lambda * zToLens * (fx.^2 + fy.^2)));
        
        % Apply thin lens phase
        focalLength = 20e-2;  % Focal length [m]
        lensPhase = exp(-1i * k0 * r.^2 / (2 * focalLength));
        uAfterLens = uAtLens .* lensPhase;
        
        % Propagate from lens to detection plane (z = f + 1 cm)
        zDetection = 21e-2;  % [m]
        uDetection = ifft2(fft2(uAfterLens) .* ...
            exp(-1i * pi * lambda * zDetection * (fx.^2 + fy.^2)));
        
        % Store intensity in cell array
        intensityCell{m1Idx} = abs(uDetection).^2;
    
    end

    % Normalize each intensity by the peak of the reference case (m1 = m2)
    fprintf('  Normalizing by reference case (m1 = %d)\n', m2);
    refPeak = max(intensityCell{1}(:));
    intensityNorm = cell(1, length(m1Values));
    for k = 1:length(intensityCell)
        intensityNorm{k} = intensityCell{k} / refPeak;
    end

    % Create full-field intensity figure
    figWhole = figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.3]);

    for m1Idx = 1:length(m1Values)
        m1 = m1Values(m1Idx);
        subplot(1, 3, m1Idx);
        imagesc(x0 * 1000, y0 * 1000, intensityNorm{m1Idx});
        axis image;
        set(gca, 'YDir', 'normal', 'FontName', 'Times New Roman', 'FontSize', 10);
        colormap hot;
        xlim([-xLim, xLim]); ylim([-yLim, yLim]);
        title(sprintf('$m_1 = %d$', m1), 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 10);
        clim([0, 1]);
        xticks(-xLim:xLim:xLim);
        yticks(-yLim:yLim:yLim);
        cb = colorbar;
        ylabel(cb, 'Normalized Intensity', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 8);
    end

    sgtitle(sprintf('Full Field Patterns for $m_2 = %d$', m2), ...
            'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 14);

    % Export full-field figure
    filename = fullfile(dirImage, sprintf('wholePattern_m2=%d.pdf', m2));
    exportgraphics(figWhole, filename);
    fprintf('  Saved: %s\n', filename);
    close(figWhole);
    
    % Clear intensity cell for next m2 iteration
    clear intensityCell intensityNorm;
end

fprintf('Done.\n');