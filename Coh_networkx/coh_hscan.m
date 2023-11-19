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
for freq=1:5
    subplot(1,5,freq)
    imagesc(squeeze(coh_mean_5freq(freq,:,:)));colorbar
end

%%
