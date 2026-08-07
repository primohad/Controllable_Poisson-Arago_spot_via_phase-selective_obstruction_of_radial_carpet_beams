% ----------------------------------------------------------
% Generate and visualize binary radial amplitude gratings
% for different spoke numbers m1.
%
% Output: ../image/m1_<m1>.pdf
%
% Notes: m1 = [10, 25, 50].
%        Amplitude grating only (no phase modulation).
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

%% ------- [ Physical Constants ] -------
% Define complex unit, wavelength, and wavenumber

i = complex(0, 1);
wl = 532e-9;    % Wavelength in meters (532 nm)
k = 2*pi/wl;    % Wavenumber [1/m]

%% ------- [ Spatial Coordinates ] -------
% Build Cartesian grid and corresponding polar coordinates

Nx = 2^10;          % Grid size along x
Ny = 2^10;          % Grid size along y

lx = 1.5e-2;        % Half-window size along x [m]
ly = 1.2e-2;        % Half-window size along y [m]

x0 = linspace(-lx, lx, Nx);
y0 = linspace(-ly, ly, Ny);
[x, y] = meshgrid(x0, y0);
[tet, r] = cart2pol(x, y);  % tet: azimuth [rad], r: radius [m]

%% ------- [ Frequency Coordinates ] -------
% Frequency axes (not used for visualization but kept for consistency)

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

%% ------- [ Main Loop: Generate Gratings ] -------
% Loop over spoke numbers and export binary amplitude gratings

m1Values = [10, 25, 50];

% Create output directory if it does not exist
dirImage = fullfile('..', 'image');
if ~exist(dirImage, 'dir')
    mkdir(dirImage);
end

for m1 = m1Values
    fprintf('Generating grating for m1 = %d\n', m1);
    
    % Create binary amplitude grating (no phase modulation)
    t1 = 0.5 * (1 + sign(cos(m1 * tet)));
    
    % Display grating pattern
    fig = figure('Color', 'white');
    imagesc(x0, y0, abs(t1).^2);
    
    % Apply consistent formatting and export
    FigureBeautificationXY(fig, lx/2, ly/2, "gray", t1, "abs");
    
    % Export with consistent naming convention
    fileName = sprintf('m1_%d.pdf', m1);
    exportgraphics(fig, fullfile(dirImage, fileName));
    
    close(fig);
end

%% ------- [ Local Functions ] -------

function FigureBeautificationXY(figHandle, limitX, limitY, colormapMode, field, dataMode)
    % Apply consistent axis limits, tick labels (mm), and color scaling
    % for intensity/phase visualization.
    %
    % Inputs:
    %   figHandle - figure handle
    %   limitX - half-width for x-axis limits (meters)
    %   limitY - half-width for y-axis limits (meters)
    %   colormapMode - 'hot' or 'gray'
    %   field - complex field data
    %   dataMode - 'abs' for intensity or 'ang' for phase
    
    axis image;
    colormap(figHandle, colormapMode);
    c = colorbar;
    
    % Set axis limits in meters
    xlim([-limitX, limitX]);
    ylim([-limitY, limitY]);
    
    % Set ticks in millimeters
    xticks([-limitX, 0, limitX]);
    set(gca, 'XTickLabel', [-limitX, 0, limitX] * 1000);
    yticks([-limitY, 0, limitY]);
    set(gca, 'YTickLabel', [-limitY, 0, limitY] * 1000);
    
    % Labels and formatting
    xlabel('$x$ (mm)', 'interpreter', 'latex');
    ylabel('$y$ (mm)', 'interpreter', 'latex');
    set(gca, 'defaulttextinterpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15);
    set(gca, 'YDir', 'normal');
    
    % Color scaling based on data mode
    if ~isempty(dataMode)
        switch dataMode
            case "abs"
                clim([0, 1]);  % Normalized intensity
                set(c, 'Ticks', [0, 1], 'TickLabelInterpreter', 'latex');
            case "ang"
                clim([-pi, pi]);
                set(c, 'TickLabelInterpreter', 'latex', ...
                    'Ticks', [-pi, pi], 'TickLabels', {'$-\pi$', '$\pi$'});
        end
    end
end