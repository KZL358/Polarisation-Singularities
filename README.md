# Polarisation Singularities

This is a companion repository to the paper [2602](https://arxiv.org/abs/), used to generate the figures, videos and simulations to compute the density of C-lines. 

## Requirements:
- MATLAB Parallel Computing Toolbox and a supported CUDA GPU.
- Script calls `gpuDevice()` and errors if no GPU is available when `useGPU=true`.

## How to Run
1. Open MATLAB and set the working directory to the folder containing the scripts.
2. Run one script at a time, e.g.
   ```matlab
   EM_Cline_density_fixed_volume.m

## Content

We split the content into three separate folders, each with their own readme file. 
* [Density](https://github.com/KZL358/Polarisation-Singularities/tree/main/Density) - This folder gives the codes necessary to compute the C-lines densities for Electromagnetic and Gravitational waves.
* [Figure](https://github.com/KZL358/Polarisation-Singularities/tree/main/Figure) - This folder gives the codes necessary to Plot the polarisation Singularities for Electromagnetic and Gravitational Waves. This is useful to reproduce Fig.(3) of the companion [paper](https://arxiv.org/abs/).
* [Video](https://github.com/KZL358/Polarisation-Singularities/tree/main/Video) - This folder gives the codes necessary to produce the video showing the deformation of polarisation singularities when perturbed. 

## Authors
* Kyan Louisia
* Claire Rigouzzo
* Sebastian Golat

