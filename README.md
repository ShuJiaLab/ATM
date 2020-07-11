# Reconstruction Demo for Airy-beam Tomographic Microscopy (ATM)
This software is distributed as accompanying software for the article Jian Wang, Xuanwen Hua, Changliang Guo, Wenhao Liu, and Shu Jia, "Airy-beam tomographic microscopy," Optica 7, 790-793 (2020),
## How to run
The demo package consists of functions and scripts written in MATLAB (MathWorks, Natick, MA). The code has been tested in MATLAB version R2018b. To simplify coding, we use DIPimage toolbox (http://www.diplib.org/).
To run the demo code for reconstruction of 3D structure:
1.	Set the current folder in MATLAB to be ATMdemo code\
2.	Open script ATMdemo.m, and using the experimental data Proj_data.mat. Set the value of the following parameters: r1_max, r2_max, r3_max (side length of reconstruction matrix), u_max (Projection side length), vpRatio (Voxel-to-pixel ratio), num_proj (number of projections), rnl (relative noise level), thetaZ (angle for the projection with z), and the file of the projection data.
3.	Run the code.
4.	The output is layer in z of the reconstruction 3D structure.
5.	The computation time depends on the r1_max, r2_max, r3_max (image size), u_max (Projection side length), vpRatio (Voxel-to-pixel ratio), num_proj (number of projections).
## Note
The program is modified based on the Open code in the project CSI:  Computational Science in Imaging, supported by the Danish Research Council for Technology and Production Sciences.
