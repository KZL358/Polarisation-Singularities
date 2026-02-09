# C-line Density Estimation (EM + GW) — MATLAB Scripts

This folder contains three MATLAB scripts that estimate **C-line density** (total C-line length per unit volume) in a **fixed physical cubic volume** for random-wave ensembles. C-lines are defined as the spatial curves where a complex scalar field vanishes:
- **Electromagnetism (EM):** $$C(\mathbf{x}) = \mathbf{E}\cdot\mathbf{E}$$
- **Gravitational waves (GW):** $$C(\mathbf{x})$$ built from the GW strain tensor components (see scripts)

Each script:
1. Builds a random field from a superposition of many plane waves (Gaussian random-wave approximation).
2. Evaluates $$C(\mathbf{x})$$ on a 3D Cartesian grid.
3. Extracts the triangulated surface **Re(C)=0** via `isosurface`.
4. Intersects that surface with **Im(C)=0** by finding sign changes across triangle edges to form short line segments.
5. Sums segment lengths to estimate total C-line length and divides by the box volume.
6. Repeats across multiple RNG seeds for statistics; reports both raw density and the rescaled quantity $$\lambda^2 \times \text{density}$$.

---

## Scripts

### 1) `EM_Cline_density_fixed_volume.m` (electromagnetic C-lines)
**Purpose:** Estimate C-line density for a random 3D complex electromagnetic field $$\mathbf{E}(\mathbf{x})$$, with
$
C(\mathbf{x}) = E_x^2 + E_y^2 + E_z^2.
$

**Key parameters (top of script):**
- `nWaves` — number of plane waves (larger → closer to Gaussian ensemble)
- `LboxHalf` — half-length of the cubic domain $$[ -L, L ]^3$$
- `lambda` — wavelength (sets $$k=2\pi/\lambda$$)
- `gridN` — grid points per axis (total points $$= \text{gridN}^3$$)
- `baseSeed`, `nSeeds`, `seeds` — reproducible random realisations

**Core function:** `compute_Cline_density_single_field(...)` (defined locally in the same file).

**Notes:**
- Builds full 3D `E(:,:,:,3)` in memory, so RAM use scales steeply with `gridN`.
- Uses single precision storage for the main field to reduce memory.

---

### 2) `GW_Cline_density_fixed_volume_cpu_blocked.m` (GW C-lines, CPU, z-blocked + parfor)
**Purpose:** Estimate C-line density for random gravitational-wave fields using a strain-tensor construction. The scalar
$
C(\mathbf{x}) = h_{11}^2 + h_{22}^2 + h_{33}^2 + 2(h_{12}^2 + h_{13}^2 + h_{23}^2)
$
is computed from complex strain components assembled from plane waves.

**Key features:**
- **z-slab blocking** (`zBlock`) to reduce peak memory.
- Pre-generates random plane-wave parameters and caches 1D exponentials `ex_all`, `ey_all`, `ez_all`.
- Uses `parfor` over blocks to speed up accumulation on multicore CPUs.
- Includes a lightweight console progress logger `consoleStageLogger`.

**Key parameters:**
- `nWaves`, `LboxHalf`, `lambda`, `gridN`
- `baseSeed`, `nSeeds`
- `zBlock` (inside function) controls memory/time tradeoff

**Core function:** `compute_GW_Cline_density_single_field(...)` (local function).

---

### 3) `GW_Cline_density_fixed_volume_gpu.m` (GW C-lines, GPU-accelerated accumulation)
**Purpose:** Same GW C-line density as (2), but uses **GPU acceleration** for the heavy tensor-field accumulation of `Re(C)` and `Im(C)`.

**Key features:**
- `useGPU = true` builds `ReC`/`ImC` on GPU in z-blocks (`zBlock` typically larger than CPU case).
- Transfers cached exponentials and per-wave coefficients to GPU.
- Geometry processing remains on CPU:
  - `isosurface` extraction
  - interpolation of `Im(C)`
  - face-loop segment extraction and length summation

**Requirements:**
- MATLAB Parallel Computing Toolbox and a supported CUDA GPU.
- Script calls `gpuDevice()` and errors if no GPU is available when `useGPU=true`.

**Core function:** `compute_GW_Cline_density_single_field(...)` (local function).

---

## How to Run

1. Open MATLAB and set the working directory to the folder containing the scripts.
2. Run one script at a time, e.g.
   ```matlab
   EM_Cline_density_fixed_volume.m
