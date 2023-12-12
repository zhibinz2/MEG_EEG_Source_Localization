clear
subj_files=[0:1:60];
for f=1:61
    tic
    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_source_hilbert
    load([num2str(subj_files(f)) '.mat'],'stroke_Cov', 'stroke_Pcov', ...
        'penalizationIn_op','penalizationOut_op','minDev_op', ...
        'subject_ID','chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');

    stroke_coh=cell(1,5);
    stroke_Pcoh=cell(1,5);
    for freq=1:5
        % convert to complex covariance and normalized to real value coherence
        cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
        stroke_coh{freq}=normalizeCSD(r2c(stroke_Cov{freq}));
        stroke_Pcoh{freq}=normalizeCSD(r2c(stroke_Pcov{freq}));
    end
    cd /home/zhibinz2/Documents/GitHub/archive/EEE_stroke_61_coh
    save([num2str(subject_ID) '.mat'],'stroke_coh', 'stroke_Pcoh', ...
        'penalizationIn_op','penalizationOut_op','minDev_op', ...
        'subject_ID','chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');
    toc
end