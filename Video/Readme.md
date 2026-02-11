This folder contains three files, one for Electromagnetic Waves and two for the Gravitational waves.
1. EM_video.m: It will create a video showing the C-lines and L-lines (polarisation singularities) for a configuration slowly going from 3 to 4 plane waves. It shows that these lines are stable under small perturbation.
2. GR_generation_frames.m: this file creates the frames, showing the C-lines and L-points for a configuration slowly going from 3 to 4 plane waves. It shows that these lines are stable under small perturbation. We separated it from the script simply because it takes longer to run. 
3. GR_video.m: Given the frames you just computed, it will plot the video.


You can change the initial parameters: 
1.nWaves: this is the final number of plane waves. There is not much difference beyond nWaves=4. 
2. lambda: it is the wavelength of the system.
3. gridN: This is the size of the grid. You can go as fine as you want, but it will become increasingly computationally heavy. 
4. rang(N): By changing the integer value of N you can change the random seed. 
5. nFrames: is the number of frames for your video. To get something smooth, we recommend a number roughly above 300.
