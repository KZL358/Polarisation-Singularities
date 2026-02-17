# C-line & L-line Visualisation (EM + GW) — MATLAB Scripts

This folder contains three MATLAB scripts that generate **videos and frames** visualising polarization singularities in random-wave fields.

We study the evolution of:

- **C-lines:** curves where a complex scalar field vanishes.
- **L-lines / L-points:** polarization singularities associated with linear polarization.

The scripts demonstrate the **stability of these singular structures under small perturbations**, by continuously deforming a configuration from **3 plane waves to 4 plane waves**.

Each workflow:

1. Builds a complex field from a controlled superposition of plane waves.
2. Computes the scalar quantity defining polarisation singularities. 
3. Extracts nodal structures on a 3D Cartesian grid.
4. Generates either individual frames or a full video showing their evolution.
5. Varies the wave configuration smoothly to illustrate structural stability.

---

## Scripts

### 1) `EM_video.m` (electromagnetic C-lines & L-lines video)

**Purpose:**  
Generates a video showing C-lines and L-lines for a complex electromagnetic field configuration that smoothly transitions from 3 to 4 plane waves.

**Notes:**
- Fully self-contained script (frame generation + video writing).
- Suitable for exploratory visualisation.
- Faster than the GR case.

---

### 2) `GR_generation_frames.m` (gravitational-wave frame generation)

**Purpose:**  
Generates individual frames showing C-lines and L-points for a gravitational-wave configuration that transitions from 3 to 4 plane waves.

This script is separated from the video creation because:
- Frame generation is computationally heavy.
- It allows reusing precomputed frames without recomputing the field.

**Output:**  
A sequence of saved frame images that can later be assembled into a video.

---

### 3) `GR_video.m` (gravitational-wave video assembly)

**Purpose:**  
Creates a video using the precomputed frames from `GR_generation_frames.m`.

**Notes:**
- Does not recompute the field.
- Only handles visualisation and video writing.
- Fast once frames are available.
- The L-points are less stable because of the interpolation procedure. We aim to fix that in future reports. 

---

## Adjustable Parameters

At the top of the scripts, you can modify:

### `nWaves`
Final number of plane waves in the configuration.

- The system smoothly transitions from 3 waves to `nWaves`.
- Little qualitative difference is observed beyond `nWaves = 4`.

---

### `lambda`
Wavelength of the system.

- Sets the spatial oscillation scale.
- Controls geometric density of structures.

---

### `gridN`
Number of grid points per spatial axis.

- Total grid points scale as `gridN^3`.
- Larger values increase resolution but significantly increase computational cost.
- Memory usage grows rapidly for 3D fields.

---

### `rng(N)`
Random seed control.

- Changing the integer `N` produces a different random realisation.
- Ensures reproducibility when fixed.

---

### `nFrames`
Number of frames in the video.

- Larger values give smoother transitions.
- Recommended: `nFrames >= 300` for visually smooth evolution.
- Directly impacts runtime.

---

## How to Run

1. Open MATLAB.
2. Set the working directory to this folder.
3. For EM:
   ```matlab
   EM_video.m
4. For GR :
   ```matlab
   GR_generation_frames.m
   GR_video.m

---

### Computational Notes 
- The EM case is relatively lightweight.
- The GR case is more computationally expensive due to its tensorial structure.
- High `gridN` and large `nFrames` substantially increase runtime.
- Separating frame generation from video assembly avoids unnecessary recomputation.
