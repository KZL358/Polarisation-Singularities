# C-line & L-line Figure Generation (EM + GW) — MATLAB Scripts

This folder contains two MATLAB scripts that generate a single 3D figure (one frame) visualising polarization singularities from a single random plane-wave realisation.

We visualise:

- **C-lines**: nodal curves where a complex scalar built from the field vanishes.
- **L-lines / L-points**: loci of perfectly linear polarization, defined via an associated n-vector constructed from the field.

Each script:

1. Builds a random complex field from a superposition of transverse plane waves.
2. Evaluates the relevant scalar and/or vector quantities on a 3D Cartesian grid.
3. Extracts C-lines (and L-lines or L-points) as intersections of zero-level sets.
4. Plots the structures inside a cubic domain.

---

## Scripts

### 1) `EM_fig.m` (electromagnetic C-lines & L-lines)

**Purpose:**  
Generate a single 3D plot of EM C-lines and L-lines from one random complex electric field **E(r)**.

- C-lines are defined by `C(r) = E·E = 0`.
- L-lines are defined by `n_E(r) = 0`, where `n_E = (1/2) Im(conj(E) × E)`.

**Method:**

- C-lines: extract `Re(C)=0` via `isosurface`, then intersect with `Im(C)=0`.
- L-lines: construct `n_E`, form two independent linear combinations, and intersect their zero sets using the same surface–contour procedure.

**Output:**  
A single figure showing:
- C-lines (red)
- L-lines (green)
- Cube outline for spatial reference

---

### 2) `GR_fig.m` (gravitational-wave C-lines & L-points)

**Purpose:**  
Generate a single 3D plot of GW C-lines and L-points from one random complex strain-tensor field **h(r)**.

- C-lines are defined by `C(r) = Tr(h^2) = 0`.
- L-points are defined by `n_h(r) = 0`, where `n_h` is constructed from `S = conj(h) * h.'`.

**Method:**

- C-lines: compute `Tr(h^2)` and intersect `Re(C)=0` with `Im(C)=0`.
- L-points: compute `n_h` on the grid, detect sign changes in each component per cell, refine candidate zeros using local minimisation of `|n_h|^2`, then cluster nearby solutions.

**Output:**  
A single figure showing:
- C-lines (red)
- L-points (green markers)
- Cube outline

---

## Key Parameters

Common parameters (top of each script):

- `nWaves` — number of random plane waves
- `lambda` — wavelength (`kMag = 2π / lambda`)
- `gridN` — grid points per axis (total points = `gridN^3`)
- `domain` — cubic domain size
- `seed`, `rngType` — reproducibility controls

Additional parameters in `GR_fig.m`:

- `valTolFrac` — tolerance fraction for accepting L-point candidates  
- `clusterTol` — clustering radius in units of grid spacing

---

## Computational Notes

- Runtime and memory scale strongly with `gridN` (as `gridN^3`).
- Field construction scales approximately linearly with `nWaves`.
- In `GR_fig.m`, L-point detection is typically the most expensive step due to full grid-cell scanning and local minimisation (`fminsearch`) of `|n_h|^2`.
- If `pagemtimes` is available, tensor contractions in the GW case are vectorised and significantly faster; otherwise a loop-based fallback is used.
- The detected L-points can appear slightly jittery depending on the choice of `valTolFrac` and `clusterTol`. We plan to improve the robustness and stability of this detection procedure in future revisions.


