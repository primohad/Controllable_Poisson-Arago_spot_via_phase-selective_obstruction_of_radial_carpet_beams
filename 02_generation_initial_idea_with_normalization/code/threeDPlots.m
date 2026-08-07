% ----------------------------------------------------------
% Visualize saved intensity matrices from ../data/
% Generates: (1) 3D beam evolution (XZ+YZ+XY planes),
%            (2) zoomed XY intensity patterns at selected z,
%            (3) longitudinal on-axis intensity traces
%
% Input: ../data/ (XZ/YZ/XY intensity matrices)
% Output: ../image/ (3D PNG, XY PDFs, longitudinal PDFs)
%
% Notes: Assumes data files exist from propagation script.
%        Uses helper function draw_surf_outline for 3D borders.
% ----------------------------------------------------------

%% ------- [ Environment Setup ] -------
% Clear workspace, close figures, compact display.
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

%% ------- [ Output Paths ] -------
% Define and create output directories if they do not exist

dirData  = fullfile('..', 'data');
% Make sure the data directory exists
if ~exist(dirData, 'dir')
    error('Data directory not found: %s', dirData);
end

dirImage = fullfile('..', 'image');

%% ------- [ Core Constants & Grid ] -------
% Define complex unit, grid sizes, and wavelength (SI units)

i = complex(0,1);
Nx = 2^11 + 1;  % Grid size along x
Ny = 2^11 + 1;  % Grid size along y
wl = 532e-9;    % Wavelength in meters (532 nm)

m1 = 10;
m2 = 10;

%% ------- [ Load Data: With Lens ] -------
% Load XZ, YZ, and selected XY intensity matrices for with-lens scenario

winHalf = 5e-3;          % [m] half-width for "with lens" scenario
prefixName = 'with lens';

f = 20e-2;
zPropagation = [linspace(0,f-(5e-2),2^8) linspace(f-(5e-2),f+(5e-2),2^8) f f-10e-2];
zPropagation = sort(unique(zPropagation)); % Ensure monotonic, unique sampling
limit = 2.5;
zTargets = [0 f-10e-2 f-5e-2, f, f+5e-2];

intensityXZPlaneYesLens = readmatrix('../data/m1 = 10 - m2 = 10 | with lens | XZ plane.txt');
intensityYZPlaneYesLens = readmatrix('../data/m1 = 10 - m2 = 10 | with lens | YZ plane.txt');

for ii = 1:numel(zTargets)
    filename = fullfile(dirData,sprintf('m1 = %d - m2 = %d | %s | f = %d cm | zPropagation = %.2f cm | Intensity XYplane.txt', ...
                                    m1, m2, prefixName, f*100, 100*zTargets(ii)));
    
    intensityXYPlaneYesLens{ii} = readmatrix(filename);
end

%% ------- [ 3D Visualization: With Lens ] -------
% Render XZ/YZ maps and embed selected XY planes as textured slices
% Top panel: 3D evolution with lens

s1 = subplot(2,1,1); hold on; colormap hot;
x0 = linspace(-winHalf, winHalf, Nx);
y0 = linspace(-winHalf, winHalf, Ny);

% XZ Plane Plot
ySurf = x0;
xSurf = zPropagation;
[XSurf,YSurf] = meshgrid(xSurf,ySurf);
ZSurf = zeros(size(YSurf));

surf1 = surf(XSurf * 100, YSurf * 1000, ZSurf - winHalf * 1000, (intensityXZPlaneYesLens), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
draw_surf_outline(surf1, 'w', 0.5);

% YZ Plane Plot
zSurf = y0;
[XSurf,ZSurf] = meshgrid(xSurf,zSurf);
YSurf = zeros(size(XSurf));

surf2 = surf(XSurf * 100, YSurf - winHalf * 1000, ZSurf * 1000, (intensityYZPlaneYesLens), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
draw_surf_outline(surf2, 'w', 0.5);

ySurf = x0;
zSurf = y0;
[YSurf,ZSurf] = meshgrid(x0,y0);
XSurf = zeros(size(YSurf));
for ii = 1:length(zTargets)
    s(ii) = surf(XSurf + zTargets(ii)*100, YSurf * 1000, ZSurf * 1000, (intensityXYPlaneYesLens{ii}), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
    draw_surf_outline(s(ii), 'w', 0.5);
end

pbaspect([5 1 1])
set(gca,'YDir','reverse')
view(-20 ,30);
set(gca,'ZTick',get(gca,'YTick'))

%% ------- [ Load Data: Without Lens ] -------
% Load XZ, YZ, and selected XY intensity matrices for free-space scenario

winHalf = 10e-3;         % [m] half-width for "without lens" scenario
prefixName = 'without lens';

zPropagation = unique(sort([linspace(0,2.5,(2^9-1)) 1.5 2.0 1.0]));
limit = winHalf * 1000;
zTargets = [0 1.0, 1.5, 2.0 2.5];

intensityXZPlaneNoLens = readmatrix('../data/m1 = 10 - m2 = 10 | without lens | XZ plane.txt');
intensityYZPlaneNoLens = readmatrix('../data/m1 = 10 - m2 = 10 | without lens | YZ plane.txt');

for ii = 1:numel(zTargets)
    filename = fullfile(dirData,sprintf('m1 = %d - m2 = %d | %s | zPropagation = %.2f cm | Intensity XYplane.txt', ...
                                    m1, m2, prefixName, 100*zTargets(ii)));
    
    intensityXYPlaneNoLens{ii} = readmatrix(filename);
end

%% ------- [ 3D Visualization: Without Lens ] -------
% Render XZ/YZ maps and embed selected XY planes as textured slices
% Bottom panel: 3D evolution without lens

s2 = subplot(2,1,2); hold on; colormap hot;
x0 = linspace(-winHalf, winHalf, Nx);
y0 = linspace(-winHalf, winHalf, Ny);

% XZ Plane Plot
ySurf = x0;
xSurf = zPropagation;
[XSurf,YSurf] = meshgrid(xSurf,ySurf);
ZSurf = zeros(size(YSurf));

surf1 = surf(XSurf * 100, YSurf * 1000, ZSurf - winHalf * 1000, (intensityXZPlaneNoLens), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
draw_surf_outline(surf1, 'w', 0.5);

% YZ Plane Plot
zSurf = y0;
[XSurf,ZSurf] = meshgrid(xSurf,zSurf);
YSurf = zeros(size(XSurf));

surf2 = surf(XSurf * 100, YSurf - winHalf * 1000, ZSurf * 1000, (intensityYZPlaneNoLens), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
draw_surf_outline(surf2, 'w', 0.5);

ySurf = x0;
zSurf = y0;
[YSurf,ZSurf] = meshgrid(x0,y0);
XSurf = zeros(size(YSurf));
for ii = 1:length(zTargets)
    s(ii) = surf(XSurf + zTargets(ii)*100, YSurf * 1000, ZSurf * 1000, (intensityXYPlaneNoLens{ii}), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
    draw_surf_outline(s(ii), 'w', 0.5);
end

pbaspect([5 1 1])
set(gca,'YDir','reverse')
view(-20 ,30);
set(gca,'ZTick',get(gca,'YTick'))

% fileName = sprintf('m1 = %d - m2 = %d | %s | 3D Plot.pdf',m1,m2,prefixName);
% exportgraphics(gcf,fullfile(dirImage,fileName))
% 
% fileName = sprintf('m1 = %d - m2 = %d | %s | 3D Plot.png',m1,m2,prefixName);
% exportgraphics(gcf,fullfile(dirImage,fileName),'Resolution',300)
% 
% fileName = sprintf('m1 = %d - m2 = %d | %s | 3D Plot.fig',m1,m2,prefixName);
% savefig(gcf,fullfile(dirImage,fileName))


% fileName = sprintf('m1 = %d - m2 = %d | 3D Plot.pdf',m1,m2);
% exportgraphics(gcf,fullfile(dirImage,fileName),'ContentType','vector')
%
set(gcf,"Units","normalized","Position",[0 0 1 1])
fileName = sprintf('m1 = %d - m2 = %d | 3D Plot.png',m1,m2);
exportgraphics(gcf,fullfile(dirImage,fileName),'Resolution',1000)

close all

%% ------- [ Transverse Planes: Zoomed XY ] -------
% Export zoomed XY intensity patterns at representative z-planes

% With lens scenario
prefixName = 'with lens';
zTargets = [f-5e-2, f, f+5e-2];

for ii = 1:numel(zTargets)
    filename = fullfile(dirData,sprintf('m1 = %d - m2 = %d | %s | f = %d cm | zPropagation = %.2f cm | Intensity XYplane.txt', ...
                                    m1, m2, prefixName, f*100, 100*zTargets(ii)));
    intensityXYPlaneYesLens{ii} = readmatrix(filename);
    figure; imagesc(x0*1000,y0*1000,intensityXYPlaneYesLens{ii});
    colormap hot; axis image
    set(gca,'YDir','normal')
    if zTargets(ii) == f
        xlim([-.35 .35])
        ylim([-.35 .35])
        xticks(-.35:.35:.35)
    else
        xlim([-1.5 1.5])
        ylim([-1.5 1.5])
        xticks(-1.5:1.5:1.5)
    end

    set(gca,'YTick',get(gca,"XTick"))
    set(gca,'YTickLabel',get(gca,'XTickLabel'))

    exportgraphics(gca,fullfile(dirImage,...
        sprintf('m1 = %d - m2 = %d | %s | f = %d cm | zPropagation = %.2f cm | Intensity XYplane.pdf', ...
                                    m1, m2, prefixName, f*100, 100*zTargets(ii))))
end

% without lens
prefixName = 'without lens';
zTargets = [1.5, 2.0 2.5];

for ii = 1:numel(zTargets)
    filename = fullfile(dirData,sprintf('m1 = %d - m2 = %d | %s | zPropagation = %.2f cm | Intensity XYplane.txt', ...
                                    m1, m2, prefixName, 100*zTargets(ii)));
    intensityXYPlaneNoLens{ii} = readmatrix(filename);
    figure; imagesc(x0*1000,y0*1000,intensityXYPlaneNoLens{ii});
    colormap hot; axis image
    set(gca,'YDir','normal')
    if ii == 1
        xlim([-3.5 3.5])
        ylim([-3.5 3.5])
        xticks(-3.5:3.5:3.5)
    elseif ii == 2
        xlim([-3.75 3.75])
        ylim([-3.75 3.75])
        xticks(-3.75:3.75:3.75)
    else
        xlim([-4.5 4.5])
        ylim([-4.5 4.5])
        xticks(-4.5:4.5:4.5)
    end

    set(gca,'YTick',get(gca,"XTick"))
    set(gca,'YTickLabel',get(gca,'XTickLabel'))

    exportgraphics(gca,fullfile(dirImage,...
        sprintf('m1 = %d - m2 = %d | %s | zPropagation = %.2f cm | Intensity XYplane.pdf', ...
                                    m1, m2, prefixName, 100*zTargets(ii))))
end

%% ------- [ Longitudinal Traces: XZ/YZ ] -------
% Plot on-axis intensity traces for XZ and YZ planes

% With lens scenario
prefixName = 'with lens';
zPropagation = [linspace(0,f-(5e-2),2^8) linspace(f-(5e-2),f+(5e-2),2^8) f f-10e-2];
zPropagation = sort(unique(zPropagation)); % Ensure monotonic, unique sampling

intensityYZPlaneYesLens = readmatrix('../data/m1 = 10 - m2 = 10 | with lens | YZ plane.txt');
intensityXZPlaneYesLens = readmatrix('../data/m1 = 10 - m2 = 10 | with lens | XZ plane.txt');

figure;
subplot(2,1,1)
plot(zPropagation * 100,intensityYZPlaneYesLens(round(Nx/2),:))
xlim([0 Inf])
ylabel('y')
xticks(0:5:25)
subplot(2,1,2)
plot(zPropagation * 100,intensityXZPlaneYesLens(round(Nx/2),:))
xlim([0 Inf])
xticks(0:5:25)
ylabel('x')

exportgraphics(gcf,fullfile(dirImage,...
    'm1 = 10 - m2 = 10 | with lens | YZ - XZ.pdf'))

% without lens
prefixName = 'without lens';
zPropagation = unique(sort([linspace(0,2.5,(2^9-1)) 1.5 2.0 1.0]));

intensityYZPlaneNoLens = readmatrix('../data/m1 = 10 - m2 = 10 | without lens | YZ plane.txt');
intensityXZPlaneNoLens = readmatrix('../data/m1 = 10 - m2 = 10 | without lens | XZ plane.txt');

figure;
subplot(2,1,1)
plot(zPropagation * 100,intensityYZPlaneNoLens(round(Nx/2),:))
xlim([0 Inf])
ylabel('y')
xticks(0:50:250)
subplot(2,1,2)
plot(zPropagation * 100,intensityXZPlaneNoLens(round(Nx/2),:))
xlim([0 Inf])
xticks(0:50:250)
ylabel('x')

exportgraphics(gcf,fullfile(dirImage,...
    'm1 = 10 - m2 = 10 | without lens | YZ - XZ.pdf'))

%% ------- [ Local Functions ] -------

function draw_surf_outline(h, color, lw)
    % Draw rectangular border around surf (outer perimeter only)
    if nargin < 2, color = 'w'; end
    if nargin < 3, lw = 1.2;   end

    X = get(h,'XData');   % size: [m n]
    Y = get(h,'YData');
    Z = get(h,'ZData');

    holdstate = ishold; hold on

    % Edges: top row, bottom row, left col, right col
    plot3(X(1,:),   Y(1,:),   Z(1,:),   '-', 'Color', color, 'LineWidth', lw);   % top edge
    plot3(X(end,:), Y(end,:), Z(end,:), '-', 'Color', color, 'LineWidth', lw);   % bottom edge
    plot3(X(:,1),   Y(:,1),   Z(:,1),   '-', 'Color', color, 'LineWidth', lw);   % left edge
    plot3(X(:,end), Y(:,end), Z(:,end), '-', 'Color', color, 'LineWidth', lw);   % right edge

    if ~holdstate, hold off, end
end