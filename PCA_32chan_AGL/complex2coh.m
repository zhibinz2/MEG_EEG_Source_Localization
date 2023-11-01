open /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util/normalizeCSD.m
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('Complex_X_op3_all.mat');

ses=1;
subj=2;
tr=11;
freq=1;
cov_c=squeeze(Complex_X_op3_all(ses,subj,tr,freq,:,:));
coh=normalizeCSD(cov_c);
imagesc(coh);colorbar;
vlim=0.001;
clim([-1*vlim vlim]);
title('coh')
subtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr) ' freq ' num2str(freq)]);

% examine the number edges
NedgeIn_coh=nan(12,2,12,5);
NedgeOut_coh=nan(12,2,12,5);
coh_all=nan(12,2,12,5,448,448);
tic
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                cov_c=squeeze(Complex_X_op3_all(ses,subj,tr,freq,:,:));
                newG=normalizeCSD(cov_c);
                coh_all(ses,subj,tr,freq,:,:)=newG;
                NedgeIn_coh(ses,subj,tr,freq) = sum(sum(logical(newG).*triu(SC,1)));
                NedgeOut_coh(ses,subj,tr,freq) =  sum(sum(logical(newG).*triu(~SC,1)));
            end
        end
    end
end
toc % 12S
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
save('coh_all.mat','coh_all','-v7.3');


%% sort through condtions
% boolean matrix average
NedgeIn_coh4=cell(4,5);
NedgeOut_coh4=cell(4,5);
tic
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                if ismember(tr,[1 2 3]); % uncouple
                    NedgeIn_coh4{1,freq}=[NedgeIn_coh4{1,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{1,freq}=[NedgeOut_coh4{1,freq} NedgeOut_coh(ses,subj,tr,freq)];
                elseif (ismember(tr,[4:6]) & subj==1) | (ismember(tr,[7:9]) & subj==2); % leading
                    NedgeIn_coh4{2,freq}=[NedgeIn_coh4{2,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{2,freq}=[NedgeOut_coh4{2,freq} NedgeOut_coh(ses,subj,tr,freq)];
                elseif (ismember(tr,[4:6]) & subj==2) | (ismember(tr,[7:9]) & subj==1); % following
                    NedgeIn_coh4{3,freq}=[NedgeIn_coh4{3,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{3,freq}=[NedgeOut_coh4{3,freq} NedgeOut_coh(ses,subj,tr,freq)];
                else ismember(tr,[10:12]);
                    NedgeIn_coh4{4,freq}=[NedgeIn_coh4{4,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{4,freq}=[NedgeOut_coh4{4,freq} NedgeOut_coh(ses,subj,tr,freq)];
                end
            end
        end
    end
end
toc

bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
NedgeIn_coh4mean=nan(4,5);
NedgeOut_coh4mean=nan(4,5);
for con=1:4
    for freq=1:5
        NedgeIn_coh4mean(con,freq)=mean(NedgeIn_coh4{con,freq});
        NedgeOut_coh4mean(con,freq)=mean(NedgeOut_coh4{con,freq});
    end
end

figure;
for con=1:4
    subplot(1,4,con)
    hold on;
    plot(1:5,NedgeIn_coh4mean(con,:)/4985,'g');
    plot(1:5,NedgeOut_coh4mean(con,:)/(99681-4985),'r');
    ylabel('n edges (prt)')
    xlabel('frequency band')
    title(['condition ' num2str(con)])
    subtitle(['percentage of edges in and outside SC']);
    xticks([1:5])
    xticklabels(bandlabels)
    ylim([0 0.2])
    hold off;
end

                
% correlation with granger




