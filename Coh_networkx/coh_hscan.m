clear
load('labels_positions.mat')

cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('coh_all.mat')

% put same frequency together in one of 5 cells
coh_5freq=cell(1,5); % 448 x 448 x 288 trials
for ses=1:12
    for subj=1:2
        for tr=1:12
            tic
            for freq=1:5
                coh_5freq{freq}=cat(3,coh_5freq{freq},squeeze(coh_all(ses,subj,tr,freq,:,:)));
            end
            toc
        end
    end
end
% mean
coh_mean_5freq=nan(5,448,448);
for freq=1:5
    coh_mean_5freq(freq,:,:)=mean(coh_5freq{freq},3);
end
% plot
figure
bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
for freq=1:5
    subplot(1,5,freq)
    imagesc(squeeze(coh_mean_5freq(freq,:,:)));colorbar
    title(bandlabels{freq})
end

coh_all_5freq=coh_5freq;
%%
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('Pcoh_boolean.mat')

% put same frequency together in one of 5 cells
coh_5freq=cell(1,5); % 448 x 448 x 288 trials
for ses=1:12
    for subj=1:2
        for tr=1:12
            tic
            for freq=1:5
                coh_5freq{freq}=cat(3,coh_5freq{freq},squeeze(Pcoh_boolean(ses,subj,tr,freq,:,:)));
            end
            toc
        end
    end
end
% mean
coh_mean_5freq=nan(5,448,448);
for freq=1:5
    coh_mean_5freq(freq,:,:)=mean(coh_5freq{freq},3);
end
% plot
figure
bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
clf
for freq=1:5
    subplot(1,5,freq)
    imagesc(squeeze(coh_mean_5freq(freq,:,:)));colorbar
    colormap('hotncold')
    title(bandlabels{freq})
end
sgtitle('boolean partial coherence avearged in source space')

Pcoh_boolean_5freq=coh_5freq;

coh_5freq=Pcoh_boolean_5freq;

% convert to zscore
coh_5all=cat(3,coh_5freq{1},coh_5freq{2},coh_5freq{3},coh_5freq{4},coh_5freq{5});
coh_5zscore=zscore(coh_5all,[],3);
% mean
coh32_mean_5freq=nan(5,448,448);
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
sgtitle('zscore: boolean partial coherence in source space averaged')

%%
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
clear
load('Pcoh_all.mat')

% put same frequency together in one of 5 cells
coh_5freq=cell(1,5); % 448 x 448 x 288 trials
for ses=1:12
    for subj=1:2
        for tr=1:12
            tic
            for freq=1:5
                coh_5freq{freq}=cat(3,coh_5freq{freq},squeeze(Pcoh_all(ses,subj,tr,freq,:,:)));
            end
            toc
        end
    end
end
% mean
coh_mean_5freq=nan(5,448,448);
for freq=1:5
    coh_mean_5freq(freq,:,:)=mean(coh_5freq{freq},3);
end
% plot
figure
bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
clf
for freq=1:5
    subplot(1,5,freq)
    imagesc(squeeze(coh_mean_5freq(freq,:,:)));colorbar
    colormap('jet')
    title(bandlabels{freq})
    clim([0 0.001])
end
sgtitle('partial coherence avearged in source space')

Pcoh_all_5freq=coh_5freq;
%%
% cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
% clear
% load('coh_all.mat')
% 
% % put same frequency together in one of 5 cells
% coh_5freq=cell(1,5); % 448 x 448 x 288 trials
% for ses=1:12
%     for subj=1:2
%         for tr=1:12
%             tic
%             for freq=1:5
%                 coh_5freq{freq}=cat(3,coh_5freq{freq},squeeze(coh_all(ses,subj,tr,freq,:,:)));
%             end
%             toc
%         end
%     end
% end
% % mean
% coh_mean_5freq=nan(5,448,448);
% for freq=1:5
%     coh_mean_5freq(freq,:,:)=mean(coh_5freq{freq},3);
% end
% % plot
% figure
% bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
% clf
% for freq=1:5
%     subplot(1,5,freq)
%     imagesc(squeeze(coh_mean_5freq(freq,:,:)));colorbar
%     colormap('jet')
%     title(bandlabels{freq})
% end
% sgtitle('ordinary coherence avearged in source space')

%%
load('/ssd/zhibin/Cleaned_sourcedata/cortical_source_data/python_lasso/Pcoh_lasso.mat')
% put same frequency together in one of 5 cells
coh_5freq=cell(1,5); % 448 x 448 x 288 trials
for ses=1:12
    for subj=1:2
        for tr=1:12
            tic
            for freq=1:5
                coh_5freq{freq}=cat(3,coh_5freq{freq},squeeze(Pcoh_lasso(ses,subj,tr,freq,:,:)));
            end
            toc
        end
    end
end
% mean
coh_mean_5freq=nan(5,448,448);
for freq=1:5
    coh_mean_5freq(freq,:,:)=mean(coh_5freq{freq},3);
end
% plot
figure
bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
clf
for freq=1:5
    subplot(1,5,freq)
    imagesc(squeeze(coh_mean_5freq(freq,:,:)));colorbar
    colormap('jet')
    title(bandlabels{freq})
    clim([0 0.001])
end
sgtitle('partial coherence by lasso averaged in source space')

Pcoh_lasso_5freq=coh_5freq;

%%
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
save('source_5freq.mat','Pcoh_lasso_5freq','coh_all_5freq','Pcoh_boolean_5freq','Pcoh_all_5freq','-v7.3');