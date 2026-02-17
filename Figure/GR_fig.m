%% GW C-lines and L-points from the same random h-field
% C-lines:  C(r) = Tr(h^2) = 0
% L-points: n_h(r) = 0, where n_h built from S = conj(h) * h.' (voxelwise)
%
% This is the GW analogue of the clean EM code structure.

clearvars; clc; close all;

%% Parameters
params.nWaves    = 10;      % number of random plane waves
params.lambda    = 1;
params.kMag      = 2*pi/params.lambda;

params.gridN     = 60;
params.domain    = 1;

params.seed      = 10;
params.rngType   = 'twister';

% L-point finder controls
params.valTolFrac = 1e-4;   % tolerance = valTolFrac * rms(|n_h|)
params.clusterTol = 1.0;    % clustering in units of grid spacing

rng(params.seed, params.rngType);

%% Grid
[x,y,z,X,Y,Z,kX,kY,kZ] = build_grid(params.gridN, -params.domain/2, params.domain/2, params.kMag);

%% Random wave ensemble (directions + amplitudes)
[K,aR,aL,phi] = generate_random_waves(params.nWaves);

%% Build h-field
tic;
h = build_random_gw_field(kX,kY,kZ, K,aR,aL,phi);
tBuildH = toc;
fprintf('h(r) built in %.3f seconds.\n', tBuildH);

%% C-lines
tic;
segmentsC = compute_gw_C_lines(X,Y,Z, x,y,z, h);
tC = toc;
fprintf('GW C-lines found in %.3f seconds.\n', tC);

%% L-points
tic;
Lpoints = compute_gw_L_points(x,y,z, h, params.valTolFrac, params.clusterTol);
tL = toc;
fprintf('GW L-points found in %.3f seconds.\n', tL);

%% Plot
plot_cube_lines_and_points(-params.domain/2, params.domain/2, segmentsC, Lpoints);

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

function [K,aR,aL,phi] = generate_random_waves(nWaves)
    % Random propagation directions (unit vectors), complex circular amplitudes,
    % plus a random global phase per wave (helps avoid accidental symmetries).
    K = randn(nWaves,3);
    K = K ./ vecnorm(K,2,2);

    aR  = (randn(nWaves,1) + 1i*randn(nWaves,1)) / sqrt(2*nWaves);
    aL  = (randn(nWaves,1) + 1i*randn(nWaves,1)) / sqrt(2*nWaves);
    phi = 2*pi*rand(nWaves,1);
end

function h = build_random_gw_field(kX,kY,kZ, K,aR,aL,phi)
    % Build complex 3x3 GW tensor field from random transverse plane waves.
    %
    % Output:
    %   h is [Nx,Ny,Nz,3,3] complex

    sz = size(kX);
    h  = complex(zeros([sz 3 3]));

    nWaves = size(K,1);

    for iWave = 1:nWaves
        khat = K(iWave,:);              % 1x3 unit direction
        [eth, eph] = transverse_basis(khat);   % SAME DEFINITION STYLE AS EM CODE

        % Plus / cross tensors (normalized like your messy code)
        Ep = (eth.'*eth - eph.'*eph)/sqrt(2);  % 3x3
        Ex = (eph.'*eth + eth.'*eph)/sqrt(2);  % 3x3

        % Circular GW helicity tensors
        ER = (Ep + 1i*Ex)/sqrt(2);
        EL = (Ep - 1i*Ex)/sqrt(2);

        H = aR(iWave)*ER + aL(iWave)*EL;       % 3x3 complex amplitude tensor

        phase = khat(1)*kX + khat(2)*kY + khat(3)*kZ + phi(iWave);
        expPhase = exp(1i*phase);              % [Nx,Ny,Nz]

        % Accumulate h(:,:,:,a,b) += H(a,b) * exp(i phase)
        for a = 1:3
            for b = 1:3
                h(:,:,:,a,b) = h(:,:,:,a,b) + H(a,b) * expPhase;
            end
        end
    end
end

function segmentsC = compute_gw_C_lines(X,Y,Z, x,y,z, h)
    % C(r) = Tr(h^2), then intersection of Re(C)=0 and Im(C)=0

    C = compute_trace_h2(h);     % complex [Nx,Ny,Nz]
    ReC = real(C);
    ImC = imag(C);

    segmentsC = intersection_curves_zeros(X,Y,Z, ReC, ImC, x,y,z);
end

function C = compute_trace_h2(h)
    % Vectorized Tr(h^2) voxelwise using pagemtimes if available.
    gridN = size(h,1);
    usePageMt = exist('pagemtimes','file') == 2;

    if usePageMt
        hReshaped  = reshape(h, [], 3, 3);                 % [N,3,3]
        h2Reshaped = pagemtimes(hReshaped, hReshaped);     % [N,3,3]
        Cvec = squeeze(h2Reshaped(:,1,1) + h2Reshaped(:,2,2) + h2Reshaped(:,3,3));
        C = reshape(Cvec, gridN, gridN, gridN);
    else
        % Fallback loops (slower, but honest)
        C = complex(zeros(gridN,gridN,gridN));
        for i = 1:gridN
            for j = 1:gridN
                for k = 1:gridN
                    Hloc = squeeze(h(i,j,k,:,:));
                    C(i,j,k) = trace(Hloc*Hloc);
                end
            end
        end
    end
end

function Lpoints = compute_gw_L_points(x,y,z, h, valTolFrac, clusterTolInGrid)
    % Compute n_h and find its zeros as L-points.

    [Nx,Ny,Nz] = compute_nh_from_h(h);

    % Interpolants
    FNx = griddedInterpolant({x, y, z}, Nx, 'linear', 'nearest');
    FNy = griddedInterpolant({x, y, z}, Ny, 'linear', 'nearest');
    FNz = griddedInterpolant({x, y, z}, Nz, 'linear', 'nearest');

    % Tolerance based on RMS magnitude
    Nh2 = Nx.^2 + Ny.^2 + Nz.^2;
    Nh_rms = sqrt(mean(Nh2(:)));
    valTol = valTolFrac * Nh_rms;

    gridN = numel(x);
    Lpoints = [];

    % fprintf('Scanning grid cells for sign changes (L-point candidates)...\n');

    for ix = 1:gridN-1
        for iy = 1:gridN-1
            for iz = 1:gridN-1

                % 8 corners of cell
                cx = [ix ix+1 ix   ix+1 ix ix+1 ix   ix+1];
                cy = [iy iy   iy+1 iy+1 iy iy   iy+1 iy+1];
                cz = [iz iz   iz   iz   iz+1 iz+1 iz+1 iz+1];

                idx = sub2ind([gridN,gridN,gridN], cx, cy, cz);

                NxC = Nx(idx); NyC = Ny(idx); NzC = Nz(idx);

                % Require sign change in each component
                if ~(min(NxC) < 0 && max(NxC) > 0), continue; end
                if ~(min(NyC) < 0 && max(NyC) > 0), continue; end
                if ~(min(NzC) < 0 && max(NzC) > 0), continue; end

                % Initial guess: cell centre
                x0 = mean(x([ix ix+1]));
                y0 = mean(y([iy iy+1]));
                z0 = mean(z([iz iz+1]));
                r0 = [x0, y0, z0];

                % Minimize |n_h|^2
                fNh2 = @(r) (FNx(r(1),r(2),r(3))).^2 + ...
                            (FNy(r(1),r(2),r(3))).^2 + ...
                            (FNz(r(1),r(2),r(3))).^2;

                opts = optimset('Display','off');
                [rStar, fval] = fminsearch(fNh2, r0, opts);

                if fval < valTol^2
                    Lpoints = [Lpoints; rStar]; %#ok<AGROW>
                end
            end
        end
    end

    % fprintf('Raw L-point candidates found: %d\n', size(Lpoints,1));

    % Deduplicate/clustering
    if ~isempty(Lpoints)
        gridSpacing = (x(end)-x(1)) / (gridN - 1);
        tolDist = clusterTolInGrid * gridSpacing;

        keep = true(size(Lpoints,1),1);
        for i = 1:size(Lpoints,1)
            if ~keep(i), continue; end
            d2 = sum((Lpoints(i+1:end,:) - Lpoints(i,:)).^2, 2);
            closeIdx = find(d2 < tolDist^2);
            keep(i + closeIdx) = false;
        end
        Lpoints = Lpoints(keep,:);
    end

    % fprintf('Distinct L-points after clustering: %d\n', size(Lpoints,1));
end

function [Nx,Ny,Nz] = compute_nh_from_h(h)
    % Build S = conj(h) * h.' (voxelwise) and n_h from its antisymmetric imaginary part.
    gridN = size(h,1);
    usePageMt = exist('pagemtimes','file') == 2;

    if usePageMt
        hReshaped = reshape(h, [], 3, 3);            % [N,3,3]
        conjH     = conj(hReshaped);                 % [N,3,3]
        hT        = permute(hReshaped, [1 3 2]);     % transpose each 3x3

        S = pagemtimes(conjH, hT);                   % [N,3,3]

        n1 = 0.5 * imag(S(:,2,3) - S(:,3,2));
        n2 = 0.5 * imag(S(:,3,1) - S(:,1,3));
        n3 = 0.5 * imag(S(:,1,2) - S(:,2,1));

        Nx = reshape(real(n1), gridN, gridN, gridN);
        Ny = reshape(real(n2), gridN, gridN, gridN);
        Nz = reshape(real(n3), gridN, gridN, gridN);
    else
        Nx = zeros(gridN,gridN,gridN);
        Ny = zeros(gridN,gridN,gridN);
        Nz = zeros(gridN,gridN,gridN);

        for i = 1:gridN
            for j = 1:gridN
                for k = 1:gridN
                    Hloc = squeeze(h(i,j,k,:,:));   % 3x3
                    S = conj(Hloc) * Hloc.';        % 3x3
                    Nx(i,j,k) = real(0.5 * imag(S(2,3) - S(3,2)));
                    Ny(i,j,k) = real(0.5 * imag(S(3,1) - S(1,3)));
                    Nz(i,j,k) = real(0.5 * imag(S(1,2) - S(2,1)));
                end
            end
        end
    end
end

function plot_cube_lines_and_points(domainMin, domainMax, segmentsC, Lpoints)
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

    if ~isempty(Lpoints)
        scatter3(Lpoints(:,1), Lpoints(:,2), Lpoints(:,3), ...
                 60, [0 0.6 0], 'filled');
    end

    camlight headlight;
    lighting gouraud;
end

function [eth, eph] = transverse_basis(khat)
    % Same spirit as your EM basis:
    %   eph = cross(khat, zhat), eth = cross(khat, eph)
    % with a deterministic choice at the poles.

    zhat = [0 0 1];
    tol  = 1e-12;

    d = dot(khat, zhat); % cos(theta)

    if abs(d) > 1 - tol
        eph = [0 1 0];
        eth = [sign(d) 0 0];
    else
        eph = cross(khat, zhat);
        eph = eph / norm(eph);

        eth = cross(khat, eph);
        eth = eth / norm(eth);
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
