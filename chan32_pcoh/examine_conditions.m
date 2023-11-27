clear
load('chan_cov_all.mat')
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



load('chan_prec_all.mat')
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
