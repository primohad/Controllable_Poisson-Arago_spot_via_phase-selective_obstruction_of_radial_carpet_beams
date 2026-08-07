% ----------------------------------------------------------
% Generate transverse intensity patterns at observation plane
% for multiple (Va, Vp) SLM settings and grating orders m1.
%
% Workflow: Apply hybrid grating on SLM aperture ->
%           Propagate to observation plane (z = 2.3m) ->
%           Export intensity with hot and gray colormaps
%
% Output: ../Va_<Va>_Vp_<Vp>/hot_propagation_m1=<m1>.pdf
%         ../Va_<Va>_Vp_<Vp>/gray_propagation_m1=<m1>.pdf
%
% Notes: Parameter sweep for comparison with experimental images.
%        m1 = [10,25,50], Va = [0.1,0.3,0.5], Vp = [pi/4,pi/2,3*pi/4].
% ----------------------------------------------------------

%% ------- [ Environment Setup ] -------
% Clear workspace, close figures, compact display

clc;
clear; close all; format compact; format bank;

%% ------- [ Physical Constants ] -------
% Define complex unit, wavelength, and wavenumber

i = complex(0, 1);
wl = 532e-9;    % Wavelength in meters (532 nm)
k = 2*pi/wl;    % Wavenumber [1/m]

%% ------- [ Spatial Coordinates ] -------
% Build Cartesian grid and corresponding polar coordinates

Nx = 2^12;          % Grid size along x
Ny = 2^12;          % Grid size along y

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

wx = 1.422e-2;      % SLM width [m]
wy = 1.066e-2;      % SLM height [m]

uPlane = zeros(size(x));
maskSLM = (x > -wx/2 & x < wx/2) & (y > -wy/2 & y < wy/2);
uPlane(maskSLM) = 1;  % Unit amplitude over SLM aperture

%% ------- [ Modulation Parameters ] -------
% Sweep parameters for amplitude visibility and phase modulation

Va = [0.1, 0.3, 0.5];           % Amplitude visibility
Vp = [pi/4, pi/2, 3*pi/4];      % Phase modulation amplitude [rad]

f = figure('Color', 'white', 'Visible', 'off');

%% ------- [ Main Loop: Sweep m1, Va, Vp ] -------
% Loop over grating orders and modulation parameters

for m1 = [10 25 50]
    fprintf('Processing m1 = %d\n', m1);

    for ii = 1:length(Va)

        for jj = 1:length(Vp)
            
            fprintf('  Va = %.1f, Vp = %.2f pi\n', Va(ii), Vp(jj)/pi);

            % Create output folder for this (Va, Vp) combination
            folderName = sprintf('Va_%.1f_Vp_%.2fpi', Va(ii), Vp(jj)/pi);
            folderPath = fullfile('..', folderName);
            if ~exist(folderPath, 'dir')
                mkdir(folderPath);
            end
            
            % grating function
            % if Va(ii) == 0
            %     t1 = exp(i .* Vp(jj) .* sign(cos(m1 .* tet)));
            % else
            %     t1 = (0.5 .* (1 + Va(ii) .* sign(cos(m1 .* tet)))) .* exp(i .* Vp(jj) .* sign(cos(m1 .* tet)));
            % end
            
            % Create hybrid radial grating with amplitude and phase modulation
            t1 = (0.5 .* (1 + Va(ii) .* sign(cos(m1 .* tet)))) ...
                .* exp(i .* Vp(jj) .* sign(cos(m1 .* tet)));
            
            % Apply grating on SLM aperture
            uGratingOnSLM = t1 .* uPlane;
            
            % Propagate from SLM plane to observation plane (z = 2.3 m)
            zp = 2.3;  % Propagation distance [m]

            uPropagated = ifft2(fft2(uGratingOnSLM) ...
                .* exp(-i .* pi .* wl * zp .*(fx.^2+fy.^2)));

            % Normalize field magnitude to unity maximum
            uPropagated = uPropagated ./ sqrt(max(max(abs(uPropagated).^2)));
            
            % Display intensity at observation plane
            imagesc(x0,y0,abs(uPropagated).^2);
            
            % Export with hot colormap
            FigureBeautificationXY(f,lx/2,ly/2,"hot",uPropagated,"abs")
            % fileName = sprintf('hot_propagation_m1=%d_Va=%.1f_Vp=%.2fpi.pdf', ...
            %        m1, Va(ii), Vp(jj)/pi);
            % exportgraphics(f, fullfile(folderPath, fileName));
            
            set(gca, "XTick", [], "YTick", [], "XLabel", [], "YLabel", [])
            colorbar('off')
            fileName = sprintf('hot_propagation_m1=%d_Va=%.1f_Vp=%.2fpi.png', ...
                   m1, Va(ii), Vp(jj)/pi);
            exportgraphics(f, fullfile(folderPath, fileName), Resolution=600);

            % Export with gray colormap
            FigureBeautificationXY(f, lx/2, ly/2, "gray", uPropagated, "abs");
            % fileName = sprintf('gray_propagation_m1=%d_Va=%.1f_Vp=%.2fpi.pdf', ...
            %                    m1, Va(ii), Vp(jj)/pi);
            % exportgraphics(f, fullfile(folderPath, fileName));
            
            set(gca, "XTick", [], "YTick", [], "XLabel", [], "YLabel", [])
            colorbar('off')
            fileName = sprintf('gray_propagation_m1=%d_Va=%.1f_Vp=%.2fpi.png', ...
                               m1, Va(ii), Vp(jj)/pi);
            exportgraphics(f, fullfile(folderPath, fileName), Resolution=600);

        end
    end
end

%% ------- [ Local Functions ] -------

function FigureBeautificationXY(figname,LimitX,LimitY,choosecolor,U,varargin)
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
    colormap(figname,choosecolor);
    c = colorbar;
    % xticks([]); yticks([]);

    % Set axis limits in meters
    xlim([-LimitX LimitX])
    ylim([-LimitY LimitY])

    % Set ticks in millimeters
    xticks(-LimitX:LimitX:LimitX); set(gca,'XTickLabel',(-LimitX:LimitX:LimitX)*1000);
    yticks(-LimitY:LimitY:LimitY); set(gca,'YTickLabel',(-LimitY:LimitY:LimitY)*1000);
    
    % Labels and formatting
    xlabel('{$x$ (mm)}','interpreter','latex')
    ylabel('{$y$ (mm)}','interpreter','latex')
    set(gca,'defaulttextinterpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'YDir','normal')
    
    % Color scaling based on data mode
    if ~isempty(varargin)
        mode = varargin{1};
        switch mode
            case "abs"
                % clim([0,ceil(max(max(abs(U).^2)))])
                % set(c,'Ticks',[0,ceil(max(max(abs(U).^2)))],'TickLabelInterpreter','latex');
                clim([0,1])
                set(c,'Ticks',[0,1],'TickLabelInterpreter','latex');
            case "ang"
                clim([-pi,pi])
                set(c,'TickLabelInterpreter','latex','Ticks',[-pi,pi],'TickLabels',{'$-\pi$','$\pi$'});
        end
    end
end