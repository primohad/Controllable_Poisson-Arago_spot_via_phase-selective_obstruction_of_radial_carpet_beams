% ----------------------------------------------------------
% Simulate diffraction from two radial gratings and analyze
% intensity patterns at the detection plane.
%
% Workflow:
%   1. Generate first hybrid grating (m1) on SLM aperture
%   2. Propagate to second grating plane (z = 2.3 m)
%   3. Apply second binary radial amplitude grating (m2)
%   4. Propagate to lens plane (z = 5 cm) and apply thin lens phase
%   5. Propagate to detection plane (z = 21 cm after lens)
%   6. Display normalized intensity pattern
%
% Input: None (all parameters defined internally)
% Output: ../image/ (intensity pattern figures)
%
% Notes: Sweeps m1 = [10, 25, 50] and m2 = [5, 10, 15, 20, 25, 30, 40, 50].
%        Detection plane at z = f + 1 cm (f = 20 cm).
%        Intensity normalized to unity maximum.
% ----------------------------------------------------------

% MATLAB R2025b

%% ------- [ Environment Setup ] -------
% Clear workspace, close figures, compact display

clearvars;
close all;
clc;
format compact;

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

dirImage  = fullfile('..', 'image');

% Create directories safely
if ~exist(dirImage, 'dir'),  mkdir(dirImage);  end

%% ------- [ Spatial Coordinates ] -------
% Build Cartesian grid and corresponding polar coordinates

Nx = 5200 + 1;          % Grid size along x
Ny = 5200 + 1;          % Grid size along y

lx = 1.5e-2;            % Half-window size along x [m]
ly = 1.2e-2;            % Half-window size along y [m]

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

u1 = zeros(size(x));
wx = 1.6e-2;            % SLM width [m]
wy = 1.2e-2;            % SLM height [m]
u1(x > -wx/2 & x < wx/2 & y > -wy/2 & y < wy/2) = 1;

%% ------- [ Optical Constants ] -------
% Define complex unit, wavelength, and wavenumber

i = complex(0, 1);
lambda = 532e-9;        % Wavelength in meters (532 nm)
k = 2*pi/lambda;        % Wavenumber [1/m]

%% ------- [ Grating Parameters ] -------
% Modulation parameters for first hybrid grating

Vp = pi/2;              % Phase modulation amplitude [rad]
Va = 0.1;               % Amplitude visibility

%% ------- [ Main Simulation Loop ] -------
% Sweep over m1 and m2 values

m1_values = [5, 10, 15, 20, 25, 30, 40, 50];
m2_values = [10, 25, 50];

% Create figure for displaying intensity patterns
fig = figure('Color', 'white', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);

for m1 = m1_values
    fprintf('Processing m1 = %d\n', m1);
    
    for m2 = m2_values
        fprintf('  m2 = %d\n', m2);
        
        % Create first hybrid radial grating (amplitude + phase)
        t1 = (0.5 * (1 + Va * sign(cos(m1 * tet)))) .* ...
             exp(1i * Vp * sign(cos(m1 * tet)));
        
        % Propagate from SLM to second grating plane (z = 2.3 m)
        zToSecondGrating = 2.3;  % [m]
        uToSecondGrating = ifft2(fft2(t1 .* u1) .* ...
            exp(-1i * pi * lambda * zToSecondGrating * (fx.^2 + fy.^2)));
        
        % Apply second binary radial amplitude grating
        t2 = 0.5 * (1 + sign(cos(m2 * tet)));
        uAtSecondGrating = uToSecondGrating .* t2;
        
        % Propagate from second grating to lens plane (z = 5 cm)
        zToLens = 5e-2;  % [m]
        uToLens = ifft2(fft2(uAtSecondGrating) .* ...
            exp(-1i * pi * lambda * zToLens * (fx.^2 + fy.^2)));
        
        % Apply thin lens phase
        f = 20e-2;  % Focal length [m]
        lensPhase = exp(-1i * k * r.^2 / (2 * f));
        uAtLens = uToLens .* lensPhase;
        
        % Propagate from lens to detection plane (z = f + 1 cm)
        zDetection = 21e-2;  % [m] (f + 1 cm)
        uDetection = ifft2(fft2(uAtLens) .* ...
            exp(-1i * pi * lambda * zDetection * (fx.^2 + fy.^2)));
        
        % Compute and normalize intensity
        intensity = abs(uDetection).^2;
        intensity = intensity ./ max(intensity(:));
        
        % Display intensity pattern
        imagesc(x0 * 1000, y0 * 1000, sqrt(intensity));
        axis image;
        set(gca, 'YDir', 'normal');
        colormap hot;
        colorbar;
        
        % Set axis limits and ticks (zoom in on central region)
        xLim = 0.35;  % [mm]
        yLim = 0.35;  % [mm]
        xlim([-xLim, xLim]);
        xticks([-xLim, 0, xLim]);
        ylim([-yLim, yLim]);
        yticks([-yLim, 0, yLim]);
        
        % Labels and title
        xlabel('x (mm)', 'FontSize', 12);
        ylabel('y (mm)', 'FontSize', 12);
        title(sprintf('m_1 = %d, m_2 = %d', m1, m2), 'FontSize', 14);
        
        % Export figure
        fileName = sprintf('intensity_m1_%d_m2_%d.png', m1, m2);
        exportgraphics(fig, fullfile(dirImage, fileName), 'Resolution', 300);
        
        % Pause briefly for visualization
        pause(0.1);
    end
end

fprintf('Simulation complete. Images saved to: %s\n', dirImage);