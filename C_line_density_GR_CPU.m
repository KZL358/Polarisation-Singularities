%% Gravitational-wave C-line density in a fixed cubic volume
% This script computes the spatial density of GW C-lines, defined as curves
% where a complex scalar invariant C(x,y,z) built from the GW strain tensor
% vanishes.
%
% For each random realisation of the GW field:
%   1) A Gaussian random wave field is synthesised from many plane waves.
%   2) The scalar field C is evaluated on a 3D Cartesian grid.
%   3) The surface Re(C)=0 is extracted as a triangulated mesh.
%   4) Intersections of this surface with Im(C)=0 define C-line segments.
%   5) The total length of all segments is divided by the box volume.
%
% The calculation is repeated for multiple random seeds to obtain statistics.

clear; clc; close all;

%% -------- Simulation parameters --------
% Number of plane waves used to approximate a Gaussian random GW field.
nWaves   = 500;

% Half-size of the cubic domain.
% The physical volume is fixed as [-LboxHalf, LboxHalf]^3.
LboxHalf = 0.5;

%% -------- Physical wavelength and grid resolution --------
% Wavelength sets the wavevector magnitude k = 2π/λ.
lambda = 0.04;

% Number of grid points per Cartesian direction (total = gridN^3).
gridN  = 600;

%% -------- Random seeds --------
% Each seed corresponds to an independent random GW field realisation.
baseSeed = 4010;
nSeeds   = 10;
seeds    = baseSeed + (0:nSeeds-1);

fprintf('Computing GW C-line density for:\n');
fprintf('  lambda = %.5f\n', lambda);
fprintf('  N_grid = %d\n', gridN);
fprintf('  nWaves = %d\n', nWaves);
fprintf('  seeds  = %d..%d (%d runs)\n', seeds(1), seeds(end), nSeeds);
fprintf('  box    = [%.2f, %.2f]^3\n\n', -LboxHalf, LboxHalf);

%% -------- Loop over random realisations --------
% Store raw and rescaled C-line densities.
lineDensity_all     = nan(nSeeds,1);
lambda2_density_all = nan(nSeeds,1);

for i = 1:nSeeds
    seed = seeds(i);
    fprintf('--- Run %d/%d (seed = %d) ---\n', i, nSeeds, seed);

    % Console logger to report progress without creating figures or GUIs.
    logp = consoleStageLogger(sprintf('seed %d', seed));

    % Compute C-line density for a single GW field realisation.
    lineDensity = compute_GW_Cline_density_single_field( ...
        nWaves, lambda, gridN, LboxHalf, seed, logp);

    lambda2_density = lambda^2 * lineDensity;

    lineDensity_all(i)     = lineDensity;
    lambda2_density_all(i) = lambda2_density;

    fprintf('  GW C-line density            = %.6e\n', lineDensity);
    fprintf('  lambda^2 * GW C-line density = %.6e\n\n', lambda2_density);
end

%% -------- Summary statistics --------
fprintf('\n==================== Summary ====================\n');
fprintf('lambda = %.5f, gridN = %d, nWaves = %d, box = [%.2f, %.2f]^3\n', ...
    lambda, gridN, nWaves, -LboxHalf, LboxHalf);

fprintf('\nPer-seed results (lambda^2 * GW C-line density):\n');
fprintf('   seed        lambda^2 * density\n');
fprintf('--------------------------------------\n');
for i = 1:nSeeds
    fprintf('%7d     %.10e\n', seeds(i), lambda2_density_all(i));
end

valid = ~isnan(lambda2_density_all);
fprintf('--------------------------------------\n');
fprintf(' valid runs: %d / %d\n', sum(valid), nSeeds);

if any(valid)
    mu  = mean(lambda2_density_all(valid));
    sig = std(lambda2_density_all(valid));
    fprintf(' mean = %.10e\n', mu);
    fprintf(' std  = %.10e\n', sig);
end

%% ============================================================
%% Local function: GW C-line density for one field realisation
%% ============================================================
function lineDensity = compute_GW_Cline_density_single_field( ...
    nWaves, lambda, gridN, LboxHalf, seed, logp)

    if nargin < 5
        error('compute_GW_Cline_density_single_field: not enough input arguments.');
    end
    if nargin < 6 || isempty(logp)
        logp = @(varargin) [];
    end

    % Fix RNG so each seed produces a reproducible GW field.
    rng(seed);
    kMag = 2*pi/lambda;

    %% -------- Block decomposition along z --------
    % The grid is processed in slabs along z to reduce peak memory usage.
    zBlock  = 8;
    nBlocks = ceil(gridN / zBlock);

    % Fixed physical volume.
    domainMin = -LboxHalf;
    domainMax =  LboxHalf;
    volume    = (domainMax - domainMin)^3;

    logp('stage', 'build grid');
    x = single(linspace(domainMin, domainMax, gridN));
    y = x;
    z = x;

    % Full coordinate grid (large but required for isosurface extraction).
    [X, Y, Z] = ndgrid(x, y, z); %#ok<ASGLU>

    %% -------- Pre-generate random plane-wave parameters --------
    % All randomness is generated once to ensure consistency and efficiency.
    kHat_all = zeros(nWaves,3);
    c1_all   = complex(zeros(nWaves,1,'single'));
    c2_all   = complex(zeros(nWaves,1,'single'));
    phi_all  = zeros(nWaves,1,'single');

    for iWave = 1:nWaves
        kHat = randn(1,3);
        kHat = kHat / norm(kHat);
        kHat_all(iWave,:) = kHat;

        c1_all(iWave)  = single(randn + 1i*randn);
        c2_all(iWave)  = single(randn + 1i*randn);
        phi_all(iWave) = single(2*pi*rand);
    end

    %% -------- Precompute per-wave tensor coefficients --------
    logp('stage', 'precompute wave coefficients');

    A = single([1, 0, 0]);

    H11_all = complex(zeros(nWaves,1,'single'));
    H22_all = complex(zeros(nWaves,1,'single'));
    H33_all = complex(zeros(nWaves,1,'single'));
    H12_all = complex(zeros(nWaves,1,'single'));
    H13_all = complex(zeros(nWaves,1,'single'));
    H23_all = complex(zeros(nWaves,1,'single'));

    ex_all   = complex(zeros(nWaves, gridN));
    ey_all   = complex(zeros(nWaves, gridN));
    ez_all   = complex(zeros(nWaves, gridN));
    ephi_all = complex(zeros(nWaves, 1));

    for iWave = 1:nWaves
        kVec = single(kMag) * single(kHat_all(iWave,:));

        if norm(cross(kVec, A)) < 1e-6
            Ause = single([0, 1, 0]);
        else
            Ause = A;
        end

        e1 = cross(kVec, Ause); e1 = e1 / norm(e1);
        e2 = cross(kVec, e1);   e2 = e2 / norm(e2);

        Ep = e1.'*e1 - e2.'*e2;
        Ex = e2.'*e1 + e1.'*e2;

        H  = (c1_all(iWave)*Ep + c2_all(iWave)*Ex) / sqrt(single(nWaves));

        H11_all(iWave) = single(H(1,1));
        H22_all(iWave) = single(H(2,2));
        H33_all(iWave) = single(H(3,3));
        H12_all(iWave) = single(H(1,2));
        H13_all(iWave) = single(H(1,3));
        H23_all(iWave) = single(H(2,3));

        ephi_all(iWave) = exp(1i * double(phi_all(iWave)));
        ex_all(iWave,:) = exp(1i * (double(kVec(1)) * double(x)));
        ey_all(iWave,:) = exp(1i * (double(kVec(2)) * double(y)));
        ez_all(iWave,:) = exp(1i * (double(kVec(3)) * double(z)));
    end

    %% -------- Assemble C field block-by-block --------
    logp('stage', 'build C per z-block');

    ReC_blk = cell(nBlocks,1);
    ImC_blk = cell(nBlocks,1);

    parfor b = 1:nBlocks
        iz0 = (b-1)*zBlock + 1;
        iz1 = min(gridN, b*zBlock);
        B   = iz1 - iz0 + 1;

        h11 = complex(zeros(gridN,gridN,B,'single'));
        h22 = h11; h33 = h11;
        h12 = h11; h13 = h11; h23 = h11;

        for iWave = 1:nWaves
            expXY = (ex_all(iWave,:).' * ey_all(iWave,:)) * ephi_all(iWave);
            slab  = expXY .* reshape(ez_all(iWave,iz0:iz1),1,1,B);

            h11 = h11 + H11_all(iWave) * slab;
            h22 = h22 + H22_all(iWave) * slab;
            h33 = h33 + H33_all(iWave) * slab;
            h12 = h12 + H12_all(iWave) * slab;
            h13 = h13 + H13_all(iWave) * slab;
            h23 = h23 + H23_all(iWave) * slab;
        end

        C = h11.^2 + h22.^2 + h33.^2 + ...
            2*(h12.^2 + h13.^2 + h23.^2);

        ReC_blk{b} = real(C);
        ImC_blk{b} = imag(C);
    end

    %% -------- Stitch full 3D fields --------
    logp('stage', 'stitch ReC/ImC');

    ReC = zeros(gridN,gridN,gridN,'single');
    ImC = zeros(gridN,gridN,gridN,'single');

    for b = 1:nBlocks
        iz0 = (b-1)*zBlock + 1;
        iz1 = min(gridN, b*zBlock);
        ReC(:,:,iz0:iz1) = ReC_blk{b};
        ImC(:,:,iz0:iz1) = ImC_blk{b};
    end

    %% -------- Extract Re(C)=0 surface --------
    logp('stage', 'isosurface Re(C)=0');

    fvRe = isosurface(double(X), double(Y), double(Z), double(ReC), 0);
    if isempty(fvRe.vertices)
        lineDensity = NaN;
        return;
    end

    %% -------- Intersect with Im(C)=0 --------
    logp('stage', 'extract C-line segments');

    FIm = griddedInterpolant({double(x), double(y), double(z)}, ...
                             double(ImC), 'linear', 'none');

    V = fvRe.vertices;
    vals = FIm(V(:,1),V(:,2),V(:,3));

    faces = fvRe.faces;
    segments = nan(3*size(faces,1),3);
    segCount = 0;

    for f = 1:size(faces,1)
        idx = faces(f,:);
        p = V(idx,:);
        v = vals(idx);

        e = [];
        if v(1)*v(2)<0, e=[e; p(1,:)+v(1)/(v(1)-v(2))*(p(2,:)-p(1,:))]; end
        if v(2)*v(3)<0, e=[e; p(2,:)+v(2)/(v(2)-v(3))*(p(3,:)-p(2,:))]; end
        if v(3)*v(1)<0, e=[e; p(3,:)+v(3)/(v(3)-v(1))*(p(1,:)-p(3,:))]; end

        if size(e,1)==2
            segments(segCount+1:segCount+2,:) = e;
            segments(segCount+3,:) = [NaN NaN NaN];
            segCount = segCount + 3;
        end
    end

    segments = segments(1:segCount,:);
    if isempty(segments)
        lineDensity = NaN;
        return;
    end

    %% -------- Length and density --------
    pts = segments(~any(isnan(segments),2),:);
    pA  = pts(1:2:end,:);
    pB  = pts(2:2:end,:);
    totalLength = sum(sqrt(sum((pB-pA).^2,2)));

    lineDensity = totalLength / volume;

    logp('done','done');
end

%% ============================================================
%% Console logger for long-running stages
%% ============================================================
function logp = consoleStageLogger(prefix)
    lastLen = 0;
    lastT   = tic;
    logp = @update;

    function update(kind,varargin)
        if strcmp(kind,'count') && toc(lastT)<0.15, return; end
        lastT = tic;

        switch kind
            case 'stage'
                line = sprintf('%s | %s ...',prefix,varargin{1});
            case 'count'
                line = sprintf('%s | %s: %d/%d',prefix,varargin{:});
            case 'done'
                line = sprintf('%s | %s',prefix,varargin{1});
            otherwise
                line = sprintf('%s | ?',prefix);
        end

        fprintf([repmat('\b',1,lastLen) '%s'],line);
        lastLen = length(line);
        if strcmp(kind,'done'), fprintf('\n'); lastLen=0; end
    end
end