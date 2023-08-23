% I picked Subject DAON, EEG code EPHD21 for testing, becasue we have this subject's 
% raw EEG data shared by Cramer, and the lesion mask, and the processed EEG
% data EPHD21

% processed EEG data
cd ../../archive/STROKE/CRAMER_directory/subacute_STROKE/EPHD21

% lesion mask
cd ../../archive/STROKE/lesion_masks/30_DAON_final_strokemask_MNI_flipped_bin.nii.gz

% raw EEG file from Cramer
cd ../../../archive/EEG_StrokePatients_n61/
load('DAON_EPHD21_20151116_1200.mat')
