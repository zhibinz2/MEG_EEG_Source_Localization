clear
subj_files=[0:1:60];
for f=1:61
    tic
    % cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_source_hilbert
    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_source_hilbert_add_delta
    load([num2str(subj_files(f)) '.mat'],'stroke_Cov', 'stroke_Pcov', ...
        'penalizationIn_op','penalizationOut_op','minDev_op', ...
        'subject_ID','chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');
    stroke_coh=cell(1,6);
    stroke_Pcoh=cell(1,6);
    for freq=1:6
        % convert to complex covariance and normalized to real value coherence
        cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
        stroke_coh{freq}=normalizeCSD(r2c(stroke_Cov{freq}));
        tmp_Cov=stroke_Pcov{freq}; tmp_Cov(isnan(tmp_Cov))=1;
        stroke_Pcoh{freq}=normalizeCSD(r2c(logical(tmp_Cov)));
    end
    cd /home/zhibinz2/Documents/GitHub/archive/EEE_stroke_61_coh
    save([num2str(subject_ID) '.mat'],'stroke_coh', 'stroke_Pcoh', ...
        'penalizationIn_op','penalizationOut_op','minDev_op', ...
        'subject_ID','chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');
    toc
end

%% troubleshooting 
f=12;freq=5;
cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_source_hilbert
load([num2str(subj_files(f)) '.mat'],'stroke_Cov', 'stroke_Pcov');
stroke_coh=cell(1,5);
stroke_Pcoh=cell(1,5);
ComplexA=r2c(stroke_Pcov{freq});

A=stroke_Pcov{5};
A(isnan(A))=1;
imagesc(logical(A));colorbar
imagesc(normalizeCSD(r2c(logical(tmp_Cov))));colorbar
r2c(A)

% convert to complex covariance and normalized to real value coherence
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
imagesc(normalizeCSD(r2c(stroke_Cov{freq})));colorbar;
imagesc(normalizeCSD(r2c(stroke_Pcov{freq})));colorbar;
