clear
%%
load('chan_cov_all.mat') % 12x2x12x5x64x64 in original order
% convert back to complex and then normalize to cross spetra
% convert to complex covariance
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
Complex_cov_all=nan(12,2,12,5,32,32);
for ses=1:numSes
    for subj=1:2
        for tr=1:12
            for freq=1:5
                Complex_cov_all(ses,subj,tr,freq,:,:)=r2c(squeeze(chan_cov_all(ses,subj,tr,freq,:,:)));
            end
        end
    end
end
% normalize to real coherence
cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
coh32_all=nan(12,2,12,5,32,32);
for ses=1:numSes
    for subj=1:2
        for tr=1:12
            for freq=1:5
                coh32_all(ses,subj,tr,freq,:,:)=normalizeCSD(squeeze(Complex_cov_all(ses,subj,tr,freq,:,:)));
            end
        end
    end
end

% sort order
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
coh32_all_sorted=nan(12,2,12,5,32,32);
for ses=1:numSes
    clear conditions sortorders
    runid = num2str(seeds(ses,:));
    load(['../../Cleaned_data/clean_' runid '.mat'],'conditions');
    % sort order
    [x,sortorder]=sort(conditions);
    for subj=1:2
        for tr=1:12
            for freq=1:5
                coh32_all_sorted(ses,subj,tr,freq,:,:) = coh32_all(ses,subj,sortorder(tr),freq,:,:);
            end
        end
    end
end


%%
load('chan_prec_all.mat')  % 12x2x12x5x64x64 in original order
% convert to complex covariance
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
Complex_prec_all=nan(12,2,12,5,32,32);
for ses=1:numSes
    for subj=1:2
        for tr=1:12
            for freq=1:5
                Complex_prec_all(ses,subj,tr,freq,:,:)=r2c(squeeze(chan_prec_all(ses,subj,tr,freq,:,:)));
            end
        end
    end
end

% normalize to real coherence
cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
coh32_all=nan(12,2,12,5,32,32);
for ses=1:numSes
    for subj=1:2
        for tr=1:12
            for freq=1:5
                coh32_all(ses,subj,tr,freq,:,:)=normalizeCSD(squeeze(Complex_prec_all(ses,subj,tr,freq,:,:)));
            end
        end
    end
end
