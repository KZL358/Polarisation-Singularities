%% C-line density computation in a fixed physical volume
% This script estimates the spatial density of C-lines (curves where
% E·E = 0) for a random 3D electromagnetic field constructed from
% a superposition of plane waves. The density is averaged over several
% random realisations (seeds) to reduce statistical fluctuations.

clear; clc; close all;

%% -------- Global simulation controls --------
% Number of plane waves used to synthesise the random field.
% Larger values better approximate a Gaussian random wave ensemble.
nWaves = 100;

% Half-length of the cubic computational domain.
% The physical box is [-LboxHalf, LboxHalf]^3 and is kept fixed
% independently of wavelength or grid resolution.
LboxHalf = 0.5;

%% -------- Wavelength and spatial resolution --------
% Wavelength used to define the wavevector magnitude.
lambda = 0.04;

% Number of grid points along each Cartesian direction.
% The total grid contains gridN^3 points.
gridN  = 500;

%% -------- Random realisations --------
% Base seed used to initialise the random number generator.
baseSeed = 2015;

% Number of independent random fields to generate.
nSeeds   = 10;

% Explicit list of seeds used for reproducibility.
seeds    = baseSeed + (0:nSeeds-1);

fprintf('Computing C-line density for:\n');
fprintf('  lambda = %.5f\n', lambda);
fprintf('  N_grid = %d\n', gridN);
fprintf('  seeds  = %d..%d (%d runs)\n', seeds(1), seeds(end), nSeeds);
fprintf('  box    = [%.2f, %.2f]^3\n\n', -LboxHalf, LboxHalf);

%% -------- Run over all random seeds --------
% Arrays storing the raw C-line density and the rescaled quantity
% lambda^2 * density for each realisation.
lineDensity_all     = nan(nSeeds,1);
lambda2_density_all = nan(nSeeds,1);

for i = 1:nSeeds
    seed = seeds(i);

    fprintf('--- Run %d/%d (seed = %d) ---\n', i, nSeeds, seed);

    % Compute the total C-line length per unit volume
    % for a single random electromagnetic field.
    lineDensity = compute_Cline_density_single_field( ...
        nWaves, lambda, gridN, LboxHalf, seed);

    % Dimensionless rescaled density commonly used for comparison across wavelengths.
    lambda2_density = lambda^2 * lineDensity;

    lineDensity_all(i)     = lineDensity;
    lambda2_density_all(i) = lambda2_density;

    fprintf('  C-line density             = %.6e\n', lineDensity);
    fprintf('  lambda^2 * C-line density  = %.6e\n\n', lambda2_density);
end

%% -------- Summary of results --------
fprintf('\n==================== Summary ====================\n');
fprintf('lambda = %.5f, gridN = %d, nWaves = %d, box = [%.2f, %.2f]^3\n', ...
    lambda, gridN, nWaves, -LboxHalf, LboxHalf);

fprintf('\nPer-seed results (lambda^2 * C-line density):\n');
fprintf('   seed        lambda^2 * density\n');
fprintf('--------------------------------------\n');
for i = 1:nSeeds
    fprintf('%7d     %.10e\n', seeds(i), lambda2_density_all(i));
end

valid = ~isnan(lambda2_density_all);
fprintf('--------------------------------------\n');
fprintf(' valid runs: %d / %d\n', sum(valid), nSeeds);

%% ============================================================
%% Local function: C-line density for one random EM field
%% ============================================================
function lineDensity = compute_Cline_density_single_field( ...
    nWaves, lambda, gridN, LboxHalf, seed)

    % Fix the random number generator so that each seed corresponds
    % to a reproducible random electromagnetic field.
    rng(seed);

    % Magnitude of the wavevector for each plane wave.
    kMag = 2*pi/lambda;

    % Define a fixed cubic domain and compute its physical volume.
    domainMin = -LboxHalf;
    domainMax =  LboxHalf;
    volume    = (domainMax - domainMin)^3;

    % Construct a uniform Cartesian grid covering the domain.
    x = single(linspace(domainMin, domainMax, gridN));
    y = x;
    z = x;
    [X, Y, Z] = ndgrid(x, y, z);

    %% -------- Random electromagnetic field construction --------
    % The field is stored as a complex-valued vector field
    % E(x,y,z) with three Cartesian components.
    E = complex(zeros(gridN, gridN, gridN, 3, 'single'));

    for iWave = 1:nWaves

        % Random propagation direction uniformly distributed on the sphere.
        kHat = randn(1,3);
        kHat = kHat / norm(kHat);
        kVec = single(kMag) * single(kHat);

        % Build an orthonormal polarisation basis transverse to k.
        if abs(dot(kHat, [0 0 1])) < 0.9
            a = [0 0 1];
        else
            a = [0 1 0];
        end
        a  = single(a);
        e1 = cross(single(kHat), a);  e1 = e1 / norm(e1);
        e2 = cross(e1, single(kHat)); e2 = e2 / norm(e2);

        % Random complex polarisation vector with unit norm.
        J  = randn(2,1) + 1i*randn(2,1);
        J  = J / norm(J);
        a1 = single(J(1));
        a2 = single(J(2));

        % Plane-wave electric field amplitude.
        % The 1/sqrt(nWaves) factor ensures finite variance
        % in the large-nWaves limit.
        Ewave = (a1*e1 + a2*e2) / sqrt(single(nWaves));

        % Random global phase.
        phi      = single(2*pi*rand);
        phase    = kVec(1)*X + kVec(2)*Y + kVec(3)*Z + phi;
        expPhase = exp(1i*phase);

        % Add this plane wave to the total field.
        E = E + expPhase .* reshape(Ewave, [1 1 1 3]);
    end

    %% -------- Complex scalar field C = E · E --------
    % C-lines are defined by the simultaneous conditions:
    %   Re(E·E) = 0
    %   Im(E·E) = 0
    Ex = E(:,:,:,1);
    Ey = E(:,:,:,2);
    Ez = E(:,:,:,3);

    C   = Ex.^2 + Ey.^2 + Ez.^2;
    ReC = real(C);
    ImC = imag(C);

    %% -------- Extract Re(C) = 0 surface --------
    % This surface is a 2D manifold embedded in 3D space.
    fvRe = isosurface(double(X), double(Y), double(Z), double(ReC), 0);

    if isempty(fvRe.vertices)
        warning('No Re(C)=0 surface found (lambda = %.3f, gridN = %d).', lambda, gridN);
        lineDensity = NaN;
        return;
    end

    %% -------- Evaluate Im(C) on the Re(C)=0 surface --------
    % Interpolate Im(C) at the vertices of the Re(C)=0 isosurface.
    FIm        = griddedInterpolant({double(x), double(y), double(z)}, ...
                                    double(ImC), 'linear', 'none');
    V          = fvRe.vertices;
    ImValsOnRe = FIm(V(:,1), V(:,2), V(:,3));

    %% -------- Locate Im(C)=0 intersections on the surface --------
    % Each triangle of the surface mesh is checked for sign changes
    % of Im(C) along its edges. Zero crossings define line segments.
    faces = fvRe.faces;
    vals  = ImValsOnRe;
    NF    = size(faces,1);

    segments = nan(3*NF, 3);
    segCount = 0;

    for f = 1:NF
        idx = faces(f,:);

        p1 = V(idx(1),:); v1 = vals(idx(1));
        p2 = V(idx(2),:); v2 = vals(idx(2));
        p3 = V(idx(3),:); v3 = vals(idx(3));

        edgePts = zeros(3,3);
        nEdge   = 0;

        if v1*v2 < 0
            t = v1/(v1 - v2);
            nEdge = nEdge + 1;
            edgePts(nEdge,:) = p1 + t*(p2 - p1);
        end
        if v2*v3 < 0
            t = v2/(v2 - v3);
            nEdge = nEdge + 1;
            edgePts(nEdge,:) = p2 + t*(p3 - p2);
        end
        if v3*v1 < 0
            t = v3/(v3 - v1);
            nEdge = nEdge + 1;
            edgePts(nEdge,:) = p3 + t*(p1 - p3);
        end

        % If exactly two edges cross Im(C)=0, a line segment is formed.
        if nEdge == 2
            segments(segCount+1,:) = edgePts(1,:);
            segments(segCount+2,:) = edgePts(2,:);
            segments(segCount+3,:) = [NaN NaN NaN];
            segCount = segCount + 3;
        end
    end

    segments = segments(1:segCount,:);
    if isempty(segments)
        warning('No C-line segments found (lambda = %.3f, gridN = %d).', lambda, gridN);
        lineDensity = NaN;
        return;
    end

    %% -------- Total C-line length --------
    % Each valid pair of points defines a short line segment.
    % Summing their lengths gives the total C-line length.
    validRows = ~any(isnan(segments), 2);
    pts = segments(validRows,:);

    pA = pts(1:2:end, :);
    pB = pts(2:2:end, :);
    segLengths = sqrt(sum((pB - pA).^2, 2));

    totalLength = sum(segLengths);

    %% -------- Line density --------
    % The C-line density is defined as total line length per unit volume.
    lineDensity = totalLength / volume;

end