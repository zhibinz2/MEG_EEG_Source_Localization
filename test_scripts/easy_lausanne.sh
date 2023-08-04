pwd
export FREESURFER_HOME=$HOME/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh


easy_lausanne \
--subject_id /home/zhibinz2/mne_data/MNE-sample-data/subjects/fsaverage \
--target_volume /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/CAMCAN_inverse.nii.gz \
--target_type diffusion \
--output_dir /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization


easy_lausanne \
--subject_id /home/zhibinz2/freesurfer/freesurfer/subjects/fsaverage \
--target_volume /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/CAMCAN_inverse.nii.gz \
--target_type diffusion \
--output_dir /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization


easy_lausanne \
--subject_id /home/zhibinz2/freesurfer/tutorial_data_20190918_1558/buckner_data/tutorial_subjs/004 \
--target_volume /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/CAMCAN_inverse.nii.gz \
--target_type diffusion \
--output_dir /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization


cd /home/zhibinz2/mne_data/MNE-sample-data/subjects/
easy_lausanne \
--subject_id fsaverage \
--target_volume /home/zhibinz2/Downloads/ds003505-download/derivatives/cmp-v3.0.3/sub-01/dwi/sub-01_desc-cmp_dwi.nii.gz \
--target_type diffusion \
--output_dir /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization

