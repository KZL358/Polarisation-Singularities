%% GW C-line density using GPU acceleration for field construction
% This script computes the density of gravitational-wave C-lines, defined as
% curves where a complex scalar C(x,y,z) constructed from the GW strain tensor
% vanishes.
%
% The algorithm proceeds as follows:
%   1) Construct a random GW field as a superposition of plane waves.
%   2) Evaluate the scalar field C over a 3D Cartesian grid.
%   3) Use GPU acceleration to assemble Re(C) and Im(C) efficiently.
%   4) Transfer the result back to CPU memory.
%   5) Extract the surface Re(C)=0 and find its intersection with Im(C)=0.
%   6) Sum the lengths of all resulting curve segments and divide by volume.
%
% GPU acceleration is used only for the heavy tensor-field accumulation.
% Isosurface extraction and geometry processing remain on the CPU.

clear; clc; close all;

%% -------- Simulation parameters --------
% Number of plane waves used to approximate a Gaussian random GW field.
nWaves      = 500;

% Enable GPU acceleration for building Re(C) and Im(C).
% All geometric post-processing is performed on the CPU.
useGPU      = true;

% Half-size of the cubic computational domain.
% The physical box is fixed as [-LboxHalf, LboxHalf]^3.
LboxHalf    = 0.5;

%% -------- Physical wavelength and grid resolution --------
% Wavelength sets the wavevector magnitude k = 2π/λ.
lambda = 0.04;

% Number of grid points per spatial dimension (total gridN^3 points).
gridN  = 600;

%% -------- Random seeds --------
% Each seed defines an independent random GW field realisation.
baseSeed   = 4010;
nSeeds     = 10;
seeds      = baseSeed + (0:nSeeds-1);

fprintf('Computing GW C-line density for:\n');
fprintf('  lambda = %.5f\n', lambda);
fprintf('  N_grid = %d\n', gridN);
fprintf('  nWaves = %d\n', nWaves);
fprintf('  seeds  = %d..%d (%d runs)\n', seeds(1), seeds(end), nSeeds);
fprintf('  box    = [%.2f, %.2f]^3\n\n', -LboxHalf, LboxHalf);

%% -------- Loop over random realisations --------
% Store raw C-line densities and λ²-rescaled values.
lineDensity_all     = nan(nSeeds,1);
lambda2_density_all = nan(nSeeds,1);

for i = 1:nSeeds
    seed = seeds(i);
    fprintf('--- Run %d/%d (seed = %d) ---\n', i, nSeeds, seed);

    % In-place console logger for long GPU computations
    logp = consoleStageLogger(sprintf('seed %d', seed));

    % Compute C-line density for this GW field realisation
    lineDensity = compute_GW_Cline_density_single_field( ...
        nWaves, lambda, gridN, useGPU, LboxHalf, seed, logp);

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
%% Local function: GW C-line density for one realisation
%% ============================================================
function lineDensity = compute_GW_Cline_density_single_field( ...
    nWaves, lambda, gridN, useGPU, LboxHalf, seed, logp)

    % Basic argument validation
    if nargin < 6
        error('compute_GW_Cline_density_single_field: not enough input arguments.');
    end
    if nargin < 7 || isempty(logp)
        logp = @(varargin) [];
    end

    % Fix RNG for reproducibility
    rng(seed);
    kMag = 2*pi/lambda;

    %% -------- z-slab decomposition --------
    % The 3D grid is processed in slabs along z to reduce GPU memory usage.
    if useGPU
        zBlock = 64;   % chosen to fit comfortably in GPU memory
    else
        zBlock = 8;
    end
    nBlocks = ceil(gridN / zBlock);

    % Fixed physical volume
    domainMin = -LboxHalf;
    domainMax =  LboxHalf;
    volume    = (domainMax - domainMin)^3;

    logp('stage', 'build grid');
    x = single(linspace(domainMin, domainMax, gridN));
    y = x;
    z = x;

    % Full coordinate grid used later for isosurface extraction
    [X, Y, Z] = ndgrid(x, y, z);

    %% -------- Random plane-wave parameters --------
    % Each plane wave has a random direction, random complex amplitudes,
    % and a random global phase.
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

    %% -------- Per-wave tensor and phase caches --------
    % To minimise GPU memory usage, only 1D exponentials are cached.
    logp('stage', 'precompute wave caches (1D ex/ey/ez)');

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

        e1 = cross(kVec, Ause);  e1 = e1 / norm(e1);
        e2 = cross(kVec, e1);    e2 = e2 / norm(e2);

        Ep = e1.'*e1 - e2.'*e2;
        Ex = e2.'*e1 + e1.'*e2;

        H  = (c1_all(iWave)*Ep + c2_all(iWave)*Ex) / sqrt(single(nWaves));

        H11_all(iWave) = single(H(1,1));
        H22_all(iWave) = single(H(2,2));
        H33_all(iWave) = single(H(3,3));
        H12_all(iWave) = single(H(1,2));
        H13_all(iWave) = single(H(1,3));
        H23_all(iWave) = single(H(2,3));

        % Phase factors stored in double for consistency with exp()
        ephi_all(iWave) = exp(1i * double(phi_all(iWave)));
        ex_all(iWave,:) = exp(1i * (double(kVec(1)) * double(x)));
        ey_all(iWave,:) = exp(1i * (double(kVec(2)) * double(y)));
        ez_all(iWave,:) = exp(1i * (double(kVec(3)) * double(z)));
    end

    %% =========================
    %% GPU computation of Re(C), Im(C)
    %% =========================
    if useGPU
        logp('stage', sprintf('GPU compute ReC/ImC blocks (zBlock=%d)', zBlock));

        % Confirm GPU availability
        try
            gpuDevice();
        catch
            error('useGPU=true but no GPU device is available.');
        end

        % Transfer cached data to GPU as complex single precision
        exG   = gpuArray(complex(single(real(ex_all)),   single(imag(ex_all))));
        eyG   = gpuArray(complex(single(real(ey_all)),   single(imag(ey_all))));
        ezG   = gpuArray(complex(single(real(ez_all)),   single(imag(ez_all))));
        ephiG = gpuArray(complex(single(real(ephi_all)), single(imag(ephi_all))));

        H11G = gpuArray(H11_all); H22G = gpuArray(H22_all); H33G = gpuArray(H33_all);
        H12G = gpuArray(H12_all); H13G = gpuArray(H13_all); H23G = gpuArray(H23_all);

        ReC_blk = cell(nBlocks,1);
        ImC_blk = cell(nBlocks,1);

        % Helper to allocate complex GPU arrays
        makeCplxZeros = @(sx,sy,sz) complex( ...
            gpuArray.zeros(sx,sy,sz,'single'), ...
            gpuArray.zeros(sx,sy,sz,'single'));

        for b = 1:nBlocks
            iz0 = (b-1)*zBlock + 1;
            iz1 = min(gridN, b*zBlock);
            B   = iz1 - iz0 + 1;

            % Accumulators for the GW strain components
            h11 = makeCplxZeros(gridN,gridN,B);
            h22 = makeCplxZeros(gridN,gridN,B);
            h33 = makeCplxZeros(gridN,gridN,B);
            h12 = makeCplxZeros(gridN,gridN,B);
            h13 = makeCplxZeros(gridN,gridN,B);
            h23 = makeCplxZeros(gridN,gridN,B);

            for iWave = 1:nWaves
                exv = exG(iWave,:).';
                eyv = eyG(iWave,:);
                expXY = (exv * eyv) * ephiG(iWave);

                ezB  = ezG(iWave, iz0:iz1);
                slab = expXY .* reshape(ezB, 1, 1, B);

                h11 = h11 + H11G(iWave) * slab;
                h22 = h22 + H22G(iWave) * slab;
                h33 = h33 + H33G(iWave) * slab;
                h12 = h12 + H12G(iWave) * slab;
                h13 = h13 + H13G(iWave) * slab;
                h23 = h23 + H23G(iWave) * slab;
            end

            C = h11.^2 + h22.^2 + h33.^2 + ...
                2*(h12.^2 + h13.^2 + h23.^2);

            ReC_blk{b} = gather(real(C));
            ImC_blk{b} = gather(imag(C));

            logp('count','zBlocks',b,nBlocks);
        end

        %% -------- Assemble full 3D fields on CPU --------
        logp('stage', 'stitch ReC/ImC blocks');

        ReC = zeros(gridN,gridN,gridN,'single');
        ImC = zeros(gridN,gridN,gridN,'single');

        for b = 1:nBlocks
            iz0 = (b-1)*zBlock + 1;
            iz1 = min(gridN, b*zBlock);
            ReC(:,:,iz0:iz1) = ReC_blk{b};
            ImC(:,:,iz0:iz1) = ImC_blk{b};
        end
    else
        error('CPU branch omitted in this GPU-focused version. Set useGPU=true.');
    end

    %% -------- Re(C)=0 surface extraction --------
    logp('stage', 'isosurface Re(C)=0');

    fvRe = isosurface(double(X), double(Y), double(Z), double(ReC), 0);
    if isempty(fvRe.vertices)
        warning('No Re(C)=0 surface found (lambda = %.3f, gridN = %d).', lambda, gridN);
        logp('done', 'no surface');
        lineDensity = NaN;
        return;
    end

    %% -------- Interpolate Im(C) on the surface --------
    logp('stage', 'interpolate Im(C) on Re(C)=0 mesh');

    FIm = griddedInterpolant({double(x), double(y), double(z)}, ...
                             double(ImC), 'linear', 'none');

    V          = fvRe.vertices;
    ImValsOnRe = FIm(V(:,1), V(:,2), V(:,3));

    %% -------- Extract C-line segments --------
    logp('stage', 'extract Im(C)=0 segments (face loop)');

    faces = fvRe.faces;
    vals  = ImValsOnRe;
    NF    = size(faces,1);

    segments = nan(3*NF, 3);
    segCount = 0;

    faceUpdateEvery = max(1, floor(NF/200));

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

        if nEdge == 2
            segments(segCount+1,:) = edgePts(1,:);
            segments(segCount+2,:) = edgePts(2,:);
            segments(segCount+3,:) = [NaN NaN NaN];
            segCount = segCount + 3;
        end

        if mod(f, faceUpdateEvery) == 0 || f == NF
            logp('count', 'faces', f, NF);
        end
    end

    segments = segments(1:segCount,:);
    if isempty(segments)
        warning('No GW C-line segments found (lambda = %.3f, gridN = %d).', lambda, gridN);
        logp('done', 'no segments');
        lineDensity = NaN;
        return;
    end

    %% -------- Compute total length and density --------
    logp('stage', 'compute total length + density');

    validRows = ~any(isnan(segments), 2);
    pts = segments(validRows,:);

    if size(pts,1) < 2
        warning('Segment extraction produced too few points (lambda = %.3f, gridN = %d).', lambda, gridN);
        logp('done', 'too few points');
        lineDensity = NaN;
        return;
    end

    pA = pts(1:2:end, :);
    pB = pts(2:2:end, :);

    segLengths  = sqrt(sum((pB - pA).^2, 2));
    totalLength = sum(segLengths);

    lineDensity = totalLength / volume;

    logp('done', 'done');
end

%% ============================================================
%% Local helper: console progress logger
%% ============================================================
function logp = consoleStageLogger(prefix)
    lastLen   = 0;
    lastPrint = tic;
    logp = @update;

    function update(kind, varargin)
        if strcmp(kind,'count') && toc(lastPrint) < 0.15
            return;
        end
        lastPrint = tic;

        switch kind
            case 'stage'
                line = sprintf('%s | %s ...', prefix, varargin{1});
            case 'count'
                line = sprintf('%s | %s: %d/%d', prefix, varargin{:});
            case 'done'
                line = sprintf('%s | %s', prefix, varargin{1});
            otherwise
                line = sprintf('%s | (unknown log event)', prefix);
        end

        fprintf([repmat('\b',1,lastLen) '%s'], line);
        lastLen = length(line);

        if strcmp(kind,'done')
            fprintf('\n');
            lastLen = 0;
        end
    end
end