% Voxel-wise DTI/DKI fitting w/ preprocessed designer data and 
% b-value and b-vector maps after GNC.
% (c) Erpeng Dai, Stanford University
% Reference: Dai E, et al. Frequency-dependent diffusion kurtosis imaging in the human brain using 
% an oscillating gradient spin echo sequence and a high-performance head-only gradient. Neuroimage 2023.

close all
clear
clc
addpath(genpath('./'))
filepath='./data/';

%% DKI
tic;
bvalname=[filepath, 'bvals'];
bvecname=[filepath, 'bvecs'];
bvalname_gnc=[bvalname, '_gnc.mat'];
bvecname_gnc=[bvecname, '_gnc.mat'];
filename=[filepath, 'data.nii.gz']; 

% % w/o GNC
% savepath_dki=[filepath, 'dki_designer_nognc/']; 
% mkdir(savepath_dki);
% cmd=['./bin/designer_mri_gnc ', ...
%     ' -DTIparams -DKIparams -WMTIparams ', ...
%     ' -fslbvec ', bvecname, ' -fslbval ', bvalname, ' ', ...
%     ...' -fslbvec_gnc ', bvecname_gnc, ' -fslbval_gnc ', bvalname_gnc, ' ', ...
%     filename,' ', savepath_dki, ...
%     ];
% [status, result]=system(cmd, '-echo');

% w/ GNC
savepath_dki=[filepath, 'dki_designer_gnc/']; 
mkdir(savepath_dki);
cmd=['./bin/designer_mri_gnc ', ...
    ' -DTIparams -DKIparams -WMTIparams ', ...
    ' -fslbvec ', bvecname, ' -fslbval ', bvalname, ' ', ...
    ' -fslbvec_gnc ', bvecname_gnc, ' -fslbval_gnc ', bvalname_gnc, ' ', ...
    filename,' ', savepath_dki, ...
    ];
[status, result]=system(cmd, '-echo');

cmd=['rm -fr ', savepath_dki, 'dwi_designer*'];
[status, result]=system(cmd, '-echo');


