%% ================================
%  PART 2 — RECORD ANIMATION TO MP4
%  ================================
clear; clc;
load('GW_animation_data.mat');

xmin = domainMin; xmax = domainMax;
ymin = domainMin; ymax = domainMax;
zmin = domainMin; zmax = domainMax;

% ---- margin parameters ----
padFrac = 0.3;   % 15% margin
dx = xmax - xmin;
dy = ymax - ymin;
dz = zmax - zmin;

xlim_pad = [xmin - padFrac*dx, xmax + padFrac*dx];
ylim_pad = [ymin - padFrac*dy, ymax + padFrac*dy];
zlim_pad = [zmin - padFrac*dz, zmax + padFrac*dz];

% ---- cube geometry ----
Vcube = [
    xmin ymin zmin;
    xmax ymin zmin;
    xmax ymax zmin;
    xmin ymax zmin;
    xmin ymin zmax;
    xmax ymin zmax;
    xmax ymax zmax;
    xmin ymax zmax ];

edges = [
    1 2; 2 3; 3 4; 4 1;
    5 6; 6 7; 7 8; 8 5;
    1 5; 2 6; 3 7; 4 8 ];

front_edges = [1 2 8 7 11 9];
back_edges  = [3 4 5 6 12 10];

%% ---- figure and axes (fixed size) ----
fig = figure('Color','w','Units','pixels','Position',[100 100 1000 800]);
ax = axes(fig);
hold(ax,'on');

axis(ax,'vis3d','equal');
view(ax,45,35);
camproj(ax,'orthographic');

xlim(ax,xlim_pad);
ylim(ax,ylim_pad);
zlim(ax,zlim_pad);

% axes padding (white border)
ax.Units = 'normalized';
ax.Position = [0.08 0.08 0.84 0.84];
ax.Visible = 'off';

% ---- draw static geometry ----
drawCubeGW(Vcube, edges, front_edges, back_edges);

Cplot = plot3(NaN,NaN,NaN,'Color',[0.82 0.01 0.11],'LineWidth',2.5);
Lplot = scatter3(NaN,NaN,NaN,30,[0.25 0.46 0.02],'filled');

camlight headlight;
lighting gouraud;

%% ---- video writer ----
videoFile = 'GW_animation_1e-3_300_frames_40_waves.mp4';
vw = VideoWriter(videoFile,'MPEG-4');
vw.FrameRate = 10;
pauseTime = 0.25;

open(vw);

%% ---- animate ONCE and record ----
for f = 1:length(ampList)
    seg = Cframes{f};
    pts = Lframes{f};

    if ~isempty(seg)
        set(Cplot,'XData',seg(:,1), ...
                  'YData',seg(:,2), ...
                  'ZData',seg(:,3));
    else
        set(Cplot,'XData',NaN,'YData',NaN,'ZData',NaN);
    end

    if ~isempty(pts)
        set(Lplot,'XData',pts(:,1), ...
                  'YData',pts(:,2), ...
                  'ZData',pts(:,3));
    else
        set(Lplot,'XData',NaN,'YData',NaN,'ZData',NaN);
    end

    title(sprintf('\\alpha = %.2f',ampList(f)),'FontSize',16);

    drawnow limitrate nocallbacks

    writeVideo(vw, getframe(fig));
    pause(pauseTime);
end

%% ---- close video cleanly ----
close(vw);
fprintf('Video written: %s\n', videoFile);

