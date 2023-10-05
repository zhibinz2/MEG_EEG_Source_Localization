%% pick the raw eeg data file from a subject
clear
cd ../../Cleaned_data/
subject_ID='20220804';
load(['clean_' subject_ID '.mat'])
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_spacing_32chan_R/


%% extract chanlocs
% cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
% cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_spacing_32chan/
chanlocs=load('chaninfo.mat')
chanlocs = chanlocs.chaninfo;

% Create channel labels
ch_labels = cell(32,1);
for c = 1:32
    ch_labels{c}=chanlocs(c).labels;
end

ch_bad=[];
dubious_chans=[];
ch_dubious=[];

%% save

Fs=2000;
% select trial and subject
trial=1;
preprocessed_eeg=dataL{trial}(:,1:32)';
save('preprocessed_eeg.mat','preprocessed_eeg','Fs','ch_dubious','ch_bad','ch_labels','chanlocs','subject_ID')


% %%
% eeg_example=preprocessed_eeg';
% time_example=1:size(eeg_example,1);
% figure;
% clf;
% plotx(time_example,eeg_example(:,[13 31])); 
% title('preprocessed-eeg');
% 
% %% test covariance
% C=cov([preprocessed_eeg]');
% figure;
% imagesc(C);colorbar();colormap('jet');clim([-1 1]);
% D=diag(C);
% mean(D)