clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% cd ../../Cleaned_data/hilbert_chancov
% for ses=1:12
%     mkdir(num2str(seeds(ses,:)));
% end

%% examine all trial
load('chan_prec_all.mat')
ses=12;
subj=2;
tr=12;
freq=5;
imagesc(squeeze(chan_prec_all(ses,subj,tr,freq,:,:))); %64x64

for ses=1:12
    for subj=1:2
        for tr=1:12
            for freq=1:5
                imagesc(squeeze(chan_prec_all(ses,subj,tr,freq,:,:)));
                colorbar;
                colormap('jet')
                clim([-1 1]);
                title(['ses ' num2str(ses) ...
                    ' subj ' num2str(subj) ...
                    ' tr ' num2str(tr) ...
                    ' freq ' num2str(freq)])
                pause(0.1);
            end 
        end
    end
end

%% convert back to complex and then normalize to cross spetra
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

% examine all trial
for ses=1:12
    for subj=1:2
        for tr=1:12
            for freq=1:5
                imagesc(squeeze(coh32_all(ses,subj,tr,freq,:,:)));
                colorbar;
                colormap('jet')
                % clim([-1 1]);
                title(['ses ' num2str(ses) ...
                    ' subj ' num2str(subj) ...
                    ' tr ' num2str(tr) ...
                    ' freq ' num2str(freq)])
                pause(0.1);
            end 
        end
    end
end


% put same frequency together in one of 5 cells
coh_5freq=cell(1,5); % 448 x 448 x 288 trials
for ses=1:12
    for subj=1:2
        for tr=1:12
            tic
            for freq=1:5
                coh_5freq{freq}=cat(3,coh_5freq{freq},logical(squeeze(coh32_all(ses,subj,tr,freq,:,:))));
            end
            toc
        end
    end
end
% mean
coh32_mean_5freq=nan(5,32,32);
for freq=1:5
    coh32_mean_5freq(freq,:,:)=mean(coh_5freq{freq},3);
end
% plot
figure
clf
bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
for freq=1:5
    subplot(1,5,freq)
    imagesc(squeeze(coh32_mean_5freq(freq,:,:)));colorbar
     colormap('jet')
    title(bandlabels{freq})
end
sgtitle('boolean partial coherence in channel space averaged')

chan_pcoh_boolean_5freq=coh_5freq;

coh_5freq=chan_pcoh_boolean_5freq;


% convert to zscore
coh_5all=cat(3,coh_5freq{1},coh_5freq{2},coh_5freq{3},coh_5freq{4},coh_5freq{5});
coh_5zscore=zscore(coh_5all,[],3);
% mean
coh32_mean_5freq=nan(5,32,32);
for freq=1:5
    coh32_mean_5freq(freq,:,:)=mean(coh_5zscore(:,:,[[1:288]+288*(freq-1)]),3);
end
% plot
figure
clf
bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
for freq=1:5
    subplot(1,5,freq)
    imagesc(squeeze(coh32_mean_5freq(freq,:,:)));colorbar
    colormap('jet')
    title(bandlabels{freq})
    clim([-0.5 0.5])
end
sgtitle('zscore: boolean partial coherence in channel space averaged')



%% examine all trial
load('chan_cov_all.mat')
ses=12;
subj=2;
tr=12;
freq=5;
imagesc(squeeze(chan_cov_all(ses,subj,tr,freq,:,:))); % 64x64
colorbar;

figure
for ses=1:12
    for subj=1:2
        for tr=1:12
            for freq=1:5
                imagesc(squeeze(chan_cov_all(ses,subj,tr,freq,:,:)));
                colorbar;
                colormap('jet')
                % clim([-1 1]);
                title(['ses ' num2str(ses) ...
                    ' subj ' num2str(subj) ...
                    ' tr ' num2str(tr) ...
                    ' freq ' num2str(freq)])
                pause(0.1);
            end 
        end
    end
end

%% convert back to complex and then normalize to cross spetra
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
% save to 'chan_5freq.mat' under /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov

% examine all trial
for ses=1:12
    for subj=1:2
        for tr=1:12
            for freq=1:5
                imagesc(squeeze(coh32_all(ses,subj,tr,freq,:,:)));
                colorbar;
                colormap('jet')
                % clim([-1 1]);
                title(['ses ' num2str(ses) ...
                    ' subj ' num2str(subj) ...
                    ' tr ' num2str(tr) ...
                    ' freq ' num2str(freq)])
                pause(0.1);
            end 
        end
    end
end


% put same frequency together in one of 5 cells
coh_5freq=cell(1,5); % 448 x 448 x 288 trials
for ses=1:12
    for subj=1:2
        for tr=1:12
            tic
            for freq=1:5
                coh_5freq{freq}=cat(3,coh_5freq{freq},squeeze(coh32_all(ses,subj,tr,freq,:,:)));
            end
            toc
        end
    end
end

% mean
coh32_mean_5freq=nan(5,32,32);
for freq=1:5
    coh32_mean_5freq(freq,:,:)=mean(coh_5freq{freq},3);
end
chan_coh_5freq=coh_5freq;
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
save('chan_5freq.mat','chan_pcoh_boolean_5freq','chan_coh_5freq','-v7.3')

% convert to zscore
coh_5all=cat(3,coh_5freq{1},coh_5freq{2},coh_5freq{3},coh_5freq{4},coh_5freq{5});
coh_5zscore=zscore(coh_5all,[],3)
% mean
coh32_mean_5freq=nan(5,32,32);
for freq=1:5
    coh32_mean_5freq(freq,:,:)=mean(coh_5zscore(:,:,[[1:288]+288*(freq-1)]),3);
end

% plot
figure
clf
bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
for freq=1:5
    subplot(1,5,freq)
    imagesc(squeeze(coh32_mean_5freq(freq,:,:)));colorbar
     colormap('jet')
    title(bandlabels{freq})
end
sgtitle('ordinary coherence in channel space averaged')
