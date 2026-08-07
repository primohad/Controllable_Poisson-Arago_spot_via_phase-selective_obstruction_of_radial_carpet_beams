% ----------------------------------------------------------
% Generate intensity/phase patterns at second grating plane
% for different (m1,m2) combinations.
%
% Workflow: Build grids -> Generate first grating (m1) ->
%           Propagate to second grating (z=2.3m) ->
%           Apply amplitude grating (m2) -> Export panels
%
% Output: ../image/abs-ang patterns at lens plane m1=<m1> m2=<m2>.pdf
%
% Notes: First pass finds global max intensity for normalization.
%        Sweeps m1 = [5,10,15,20,25,30,40,50], m2 = [10,25,50].
% ----------------------------------------------------------

%% ------- [ Environment Setup ] -------
% Clear workspace, close figures, compact display

clc
clear; close all; format compact; format bank;

%% ------- [ Working Directory Sync ] -------
% Auto-change current folder to the active script's location.
disp('Auto change folder enabled.');
activeFile = matlab.desktop.editor.getActiveFilename;
if ~isempty(activeFile)
    cd(fileparts(activeFile));
    fprintf('Changed folder to: %s\n', pwd);
end

%% ------- [ Optical Constants ] -------
% Define complex unit, wavelength, and wavenumber

i = complex(0, 1);
wl = 532e-9;    % Wavelength in meters (532 nm)
k = 2*pi/wl;    % Wavenumber [1/m]

%% ------- [ Spatial Coordinates ] -------

Nx = 2^10 + 1;
Ny = 2^10 + 1;

% Spatial window size (meters)
lx = 1.5e-2;    % Half-width along x [m]
ly = 1.2e-2;    % Half-width along y [m]
x0 = linspace(-lx, lx, Nx);
y0 = linspace(-ly, ly, Ny);
[x,y] = meshgrid(x0,y0);
[tet,r] = cart2pol(x,y);

%% ------- [ Frequency Coordinates ] -------

Lx = max(x0)-min(x0); % window width
Ly = max(y0)-min(y0); % window height

dx = Lx/Nx;
dy = Ly/Ny;

dfx = 1/(Nx*dx); % fundamental frequency in frequency-domain
dfy = 1/(Ny*dy); % fundamental frequency in frequency-domain

fx0 = dfx.*(-Nx/2:Nx/2-1); % frequency (along horizontal direction)
fy0 = dfy.*(-Ny/2:Ny/2-1); % frequency (along vertical direction)

fx0 = fftshift(fx0);
fy0 = fftshift(fy0);

[fx,fy] = meshgrid(fx0,fy0); % creating fx-fy plane

%% ------- [ Input Field: SLM Aperture ] -------

% SLM aperture dimensions (meters)
wx = 1.6e-2;    % Width along x [m]
wy = 1.2e-2;    % Width along y [m]

% plane wave on SLM
uPlane = zeros(size(x));
uPlane(x > -wx/2 & x < wx/2 & y > -wy/2 & y < wy/2) = 1;

%% ------- [ First Pass: Global Max Intensity ] -------
% Propagate without visualization to find global maximum for normalization
fprintf('First pass: finding global maximum intensity...\n');

globalMaxIntensity = 0;
zPropagation = 2.3;    % Distance to second grating [m]

vp = pi/2;      % Phase modulation amplitude [rad]
va = 0.1;       % Visibility (amplitude modulation depth)

for m1 = [5,10,15,20,25,30,40,50]
    disp(['m1 = ', num2str(m1), ' (first pass)']);

    % grating function
    t1 = (0.5 .* (1 + va .* sign(cos(m1 .* tet)))) .* exp(i .* vp .* sign(cos(m1 .* tet)));
    
    % grating on SLM
    uGratingOnSLM = t1 .* uPlane;
    
    % Propagate from SLM plane to second-grating plane (z = 2.3 m)
    uPropagateBefore = ifft2(fft2(uGratingOnSLM) .* exp(-i .* pi .* wl * zPropagation .*(fx.^2+fy.^2)));
    
    % Apply binary radial amplitude grating (m2)
    for m2 = [10,25,50]
        disp(['m2 = ', num2str(m2), ' (first pass)']);
    
        t2 = 0.5 .* (1 + sign(cos(m2 * tet))); % binary amplitude grating
        uAtSecondGrating = uPropagateBefore .* t2;
        
        % Get maximum intensity for this iteration
        currentMax = max(abs(uAtSecondGrating(:)).^2);
        
        % Update global maximum
        if currentMax > globalMaxIntensity
            globalMaxIntensity = currentMax;
        end
    end
end

disp(['Global maximum intensity found: ', num2str(globalMaxIntensity)]);
disp(['Ceiling value for colorbar: ', num2str(ceil(globalMaxIntensity))]);

%% ------- [ Main Loop: Sweep m1 and m2 ] -------

vp = pi/2;      % Phase modulation amplitude [rad]
va = 0.1;       % Visibility (amplitude modulation depth)

for m1 = [5, 10, 15, 20, 25, 30, 40, 50]
    disp(['m1 = ', num2str(m1)]);

    % grating function
    t1 = (0.5 .* (1 + va .* sign(cos(m1 .* tet)))) .* exp(i .* vp .* sign(cos(m1 .* tet)));
    
    % grating on SLM
    uGratingOnSLM = t1 .* uPlane;
    
    % Propagate from SLM plane to second-grating plane (z = 2.3 m)
    zPropagation = 2.3; % m
    uPropagateBefore = ifft2(fft2(uGratingOnSLM) .* exp(-i .* pi .* wl * zPropagation .*(fx.^2+fy.^2)));
    
    % Apply binary radial amplitude grating (m2)
    for m2 = [10,25,50]
        disp(['m2 = ', num2str(m2)]);
    
        t2 = 0.5 .* (1 + sign(cos(m2 * tet))); % binary amplitude grating
        
        uAtSecondGrating = uPropagateBefore .* t2 / sqrt(globalMaxIntensity);
    
        % Create figure with intensity (top) and phase (bottom) panels
        fig = figure;

        % Intensity panel
        s(1) = subplot(2,1,1);
        imagesc(x0,y0, double(abs(uAtSecondGrating).^2));
        FigureModificationLimited(s(1), [-wy/2 wy/2], [-wy/2 wy/2], "abs", abs(uAtSecondGrating).^2)
        
        % Phase panel
        s(2) = subplot(2,1,2);
        imagesc(x0,y0,angle(uAtSecondGrating))
        FigureModificationLimited(s(2), [-wy/2 wy/2], [-wy/2 wy/2], "ang", phase)
        
        % Export intensity/phase panels
        exportgraphics(fig,['../image/' 'abs-ang patterns at second grating plane m1=' num2str(m1) ' m2=' num2str(m2) '.pdf'])
    end
end

%% ------- [ Local Functions ] -------

function FigureModificationLimited(name, x0boundary, y0boundary, mode, func)
    % Apply consistent axis limits, tick labels (mm), and color scaling
    % for intensity/phase visualization.
    %
    % Inputs:
    %   axHandle - axes handle
    %   xBoundary - [xmin, xmax] for x-axis limits (meters)
    %   yBoundary - [ymin, ymax] for y-axis limits (meters)
    %   mode - 'abs' for intensity (hot colormap) or 'ang' for phase (jet colormap)
    %   data - matrix for color scaling
    
    axis image;
    
    if mode == "abs"
        colormap(name, "hot");
        c = colorbar;

        % clim([0,1])
        % set(c,'TickLabelInterpreter','latex','Ticks',[0,1]);

        clim([floor(min(min(func))),ceil(max(max(func)))]);
        set(c,'TickLabelInterpreter','latex','Ticks',[floor(min(min(func))),ceil(max(max(func)))])
    elseif mode == "ang"
        colormap(name, "jet");
        c = colorbar;
        
        clim([-pi,pi])
        set(c,'TickLabelInterpreter','latex','Ticks',[-pi,pi],'TickLabels',{'$-\pi$','$\pi$'});
    
        % % Optional alternative phase scaling kept for experiments; default uses [-pi, pi]
        % if (floor(min(min(func))) == 0 && ceil(max(max(func))) == 0)
        %     clim([-pi,pi])
        %     set(c,'TickLabelInterpreter','latex','Ticks',[-pi,pi],'TickLabels',{'$-\pi$','$\pi$'});
        % 
        % else
        %     clim([floor(min(min(func))),ceil(max(max(func)))]);
        %     set(c,'TickLabelInterpreter','latex','Ticks',[floor(min(min(func))),ceil(max(max(func)))], ...
        %         'TickLabels',{['$-' num2str(floor(min(min(func)))/pi) '\pi$'],'$' num2str(ceil(max(max(func)))/pi) '\pi$'})
        % end
    end
    
    xlim([min(x0boundary), max(x0boundary)]);
    ylim([min(y0boundary), max(y0boundary)]);

    xticks([min(x0boundary), max(x0boundary)]); set(name,'XTickLabel',[min(x0boundary), max(x0boundary)]*1000);
    yticks([min(y0boundary), max(y0boundary)]); set(name,'YTickLabel',[min(y0boundary), max(y0boundary)]*1000);

    xlabel('{$x$(mm)}','interpreter','latex')
    ylabel('{$y$(mm)}','interpreter','latex')

    set(name,'defaulttextinterpreter','latex')
    set(name,'TickLabelInterpreter','latex','FontSize',15)
    set(name,'YDir','normal')
end