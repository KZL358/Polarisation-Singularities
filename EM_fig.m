%% EM C-lines and L-lines from the same random E-field
% C-lines:  C(r) = E·E = 0   (perfectly circular polarization)
% L-lines:  n_E(r) = 0       (perfectly linear polarization)
% where n_E = (1/2) Im(E* x E) = P x Q in Berry's notation

clearvars; clc; close all;

%% Parameters
params.nWaves    = 100;
params.lambda    = 1;
params.kMag      = 2*pi/params.lambda;

params.gridN     = 60;
params.domain    = 1;

params.seed      = 23;
params.rngType   = 'twister';

rng(params.seed, params.rngType);

%% Grid
[x,y,z,X,Y,Z,kX,kY,kZ] = build_grid(params.gridN, -params.domain/2, params.domain/2, params.kMag);

%% Random wave ensemble (pre-generated)
[K,aR,aL] = generate_random_waves(params.nWaves);

%% Build E-field
tic;
[Ex,Ey,Ez] = build_random_field(kX,kY,kZ, K,aR,aL);
tBuildE = toc;
fprintf('E(r) built in %.3f seconds.\n', tBuildE);

%% C-lines
tic;
segmentsC = compute_C_lines(X,Y,Z, x,y,z, Ex,Ey,Ez);
tC = toc;
fprintf('EM C lines found in %.3f seconds.\n', tC);

%% L-lines
tic;
segmentsL = compute_L_lines(X,Y,Z, x,y,z, Ex,Ey,Ez);
tL = toc;
fprintf('EM L lines found in %.3f seconds.\n', tL);

%% Plot
plot_cube_and_lines(-params.domain/2, params.domain/2, segmentsC, segmentsL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions

function [x,y,z,X,Y,Z,kX,kY,kZ] = build_grid(gridN, domainMin, domainMax, kMag)
    x = linspace(domainMin, domainMax, gridN);
    y = linspace(domainMin, domainMax, gridN);
    z = linspace(domainMin, domainMax, gridN);

    [X,Y,Z] = ndgrid(x,y,z);

    kX = kMag * X;
    kY = kMag * Y;
    kZ = kMag * Z;
end

function [K,aR,aL] = generate_random_waves(nWaves)
    K = randn(nWaves,3);
    K = K ./ vecnorm(K,2,2);

    aR = (randn(nWaves,1) + 1i*randn(nWaves,1)) / sqrt(2*nWaves);
    aL = (randn(nWaves,1) + 1i*randn(nWaves,1)) / sqrt(2*nWaves);
end

function [Ex,Ey,Ez] = build_random_field(kX,kY,kZ, K,aR,aL)
    % Build complex 3D vector field from random transverse plane waves.
    %
    % kX,kY,kZ are 3D arrays containing kMag*X etc.
    % K is [nWaves,3], aR,aL are [nWaves,1].

    sz = size(kX);
    Ex = complex(zeros(sz));
    Ey = complex(zeros(sz));
    Ez = complex(zeros(sz));

    nWaves = size(K,1);

    for iWave = 1:nWaves
        khat = K(iWave,:);
        [eR, eL] = transverse_basis(khat);  

        p = aR(iWave)*eR + aL(iWave)*eL;    % 1x3 complex

        phase = khat(1)*kX + khat(2)*kY + khat(3)*kZ;
        expPhase = exp(1i*phase);

        Ex = Ex + expPhase * p(1);
        Ey = Ey + expPhase * p(2);
        Ez = Ez + expPhase * p(3);
    end
end

function segmentsC = compute_C_lines(X,Y,Z, x,y,z, Ex,Ey,Ez)
    C   = Ex.^2 + Ey.^2 + Ez.^2;
    ReC = real(C);
    ImC = imag(C);
    segmentsC = intersection_curves_zeros(X,Y,Z, ReC, ImC, x,y,z);
end

function segmentsL = compute_L_lines(X,Y,Z, x,y,z, Ex,Ey,Ez)
    Nx = 0.5 * imag( conj(Ey).*Ez - conj(Ez).*Ey );
    Ny = 0.5 * imag( conj(Ez).*Ex - conj(Ex).*Ez );
    Nz = 0.5 * imag( conj(Ex).*Ey - conj(Ey).*Ex );

    % choose two independent linear combinations
    a = [1; 0; 0];
    b = [0; 1; 1] / sqrt(2);

    f = a(1)*Nx + a(2)*Ny + a(3)*Nz;
    g = b(1)*Nx + b(2)*Ny + b(3)*Nz;

    segmentsL = intersection_curves_zeros(X,Y,Z, f, g, x,y,z);
end

function plot_cube_and_lines(domainMin, domainMax, segmentsC, segmentsL)
    figure('Color','w'); hold on;

    xmin = domainMin; xmax = domainMax;
    ymin = domainMin; ymax = domainMax;
    zmin = domainMin; zmax = domainMax;

    V = [
        xmin ymin zmin;
        xmax ymin zmin;
        xmax ymax zmin;
        xmin ymax zmin;
        xmin ymin zmax;
        xmax ymin zmax;
        xmax ymax zmax;
        xmin ymax zmax
    ];

    edges = [
        1 2; 2 3; 3 4; 4 1;
        5 6; 6 7; 7 8; 8 5;
        1 5; 2 6; 3 7; 4 8
    ];

    az = 45; el = 35;
    view(az, el);
    camproj orthographic;

    axis equal
    xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]);
    set(gca,'Visible','off','XTick',[],'YTick',[],'ZTick',[]);
    box off;

    front_edges = [1 2 8 7 11 9];
    back_edges  = [3 4 5 6 12 10];

    for ei = 1:size(edges,1)
        p1 = V(edges(ei,1),:);
        p2 = V(edges(ei,2),:);

        if ismember(ei, front_edges)
            style = '-';  col = [0 0 0];       lw = 1.6;
        elseif ismember(ei, back_edges)
            style = '--'; col = [0.7 0.7 0.7]; lw = 1.0;
        else
            style = '-';  col = [0.55 0.55 0.55]; lw = 1.0;
        end

        plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], ...
              'LineStyle', style, 'Color', col, 'LineWidth', lw);
    end

    if ~isempty(segmentsC)
        plot3(segmentsC(:,1), segmentsC(:,2), segmentsC(:,3), ...
              'Color', [0.9 0 0], 'LineWidth', 4.0);
    end

    if ~isempty(segmentsL)
        plot3(segmentsL(:,1), segmentsL(:,2), segmentsL(:,3), ...
              'Color', [0 0.6 0], 'LineWidth', 4.0);
    end

    camlight headlight;
    lighting gouraud;
end

function [e1, e2] = transverse_basis(khat, mode)
    if nargin < 2, mode = 'circular'; end

    zhat = [0 0 1];

    eph = cross(khat, zhat);
    eph = eph / norm(eph);

    eth = cross(khat, eph);
    eth = eth / norm(eth);

    switch lower(mode)
        case 'linear'
            e1 = eth;
            e2 = eph;
        case 'circular'
            e1 = (eth + 1i*eph)/sqrt(2);
            e2 = (eth - 1i*eph)/sqrt(2);
        otherwise
            error('transverse_basis: unknown mode "%s"', mode);
    end
end

function segments = intersection_curves_zeros(X,Y,Z, f, g, x,y,z, varargin)
    if isempty(varargin)
        f0 = 0; g0 = 0;
    elseif numel(varargin) == 2
        f0 = varargin{1};
        g0 = varargin{2};
    else
        error('intersection_curves_zeros: expected 0 or 2 optional args (f0,g0).');
    end

    fv = isosurface(X, Y, Z, f, f0);
    if isempty(fv.vertices)
        segments = zeros(0,3);
        return
    end

    V = fv.vertices;
    F = fv.faces;

    Ginterp = griddedInterpolant({x,y,z}, g, 'linear', 'nearest');
    gv = Ginterp(V(:,1), V(:,2), V(:,3));

    segments = contour_on_trimesh(V, F, gv, g0);
end

function segments = contour_on_trimesh(V, F, vals, level)
    vals = vals(:) - level;

    goodFace = all(isfinite(vals(F)), 2);
    F = F(goodFace,:);

    eps0 = 1e-12;

    segCell = cell(size(F,1),1);
    nSeg = 0;

    edges = [1 2; 2 3; 3 1];

    for fi = 1:size(F,1)
        idx = F(fi,:);
        p = V(idx,:);
        v = vals(idx);

        v(abs(v) < eps0) = 0;

        pts = zeros(2,3);
        n = 0;

        for e = 1:3
            i = edges(e,1); j = edges(e,2);
            vi = v(i); vj = v(j);

            if vi == 0 && vj == 0
                continue
            elseif vi == 0
                n = n + 1; if n <= 2, pts(n,:) = p(i,:); end
            elseif vj == 0
                n = n + 1; if n <= 2, pts(n,:) = p(j,:); end
            elseif vi * vj < 0
                t = vi / (vi - vj);
                n = n + 1;
                if n <= 2
                    pts(n,:) = p(i,:) + t*(p(j,:) - p(i,:));
                end
            end
        end

        if n == 2
            nSeg = nSeg + 1;
            segCell{nSeg} = [pts(1,:); pts(2,:); NaN NaN NaN];
        end
    end

    if nSeg == 0
        segments = zeros(0,3);
    else
        segments = vertcat(segCell{1:nSeg});
    end
end