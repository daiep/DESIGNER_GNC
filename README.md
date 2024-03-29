# DESIGNER_GNC
Voxel-wise DKI/DTI fitting script with Gradient Non-linearity Corrected (GNC) b-value and b-vector maps, based on the DESIGNER toolbox (_Ades-Aron B, et al. NIMG 2018_).
## Getting Started
Run DESIGNER_DKI_GNC.m for voxel-wise DKI/DTI fitting, with synthesized multi-shell diffusion data (already masked and only for demo purpose).
  
Refer to the example data for the data format requirement:  
1.  b0/data/nodif_brain_mask.nii.gz, bvals, and bvecs: Preprocessed data with the recommended DESIGNER pipeline, including corrections for thermal noise, Gibbs ringing, Rician biases, susceptibility–induced geometric distortions, eddy-current-induced spatial distortions, and interslice motion.
2.  bvals_gnc.mat and bvecs_gnc.mat: Voxel-wise b-value and b-vector maps after GNC. Format: #-of-b-val (or b-vec) * 3 (for b-vec)/1 (for b-val) * RO * PE * SL. 
## Citations
If you find this script useful in your research, please cite:  

Dai E, Zhu A, Yang GK, et al. Frequency-dependent diffusion kurtosis imaging in the human brain using an oscillating gradient spin echo sequence and a high-performance head-only gradient. Neuroimage 2023; 120328.  
  
and  
### Original DESIGNER paper/toolbox
Ades-Aron B, Veraart J, Kochunov P, et al. Evaluation of the accuracy and precision of the diffusion parameter EStImation with Gibbs and NoisE removal pipeline. Neuroimage 2018 (183), 532-543. 

https://github.com/NYU-DiffusionMRI/DESIGNER.git
