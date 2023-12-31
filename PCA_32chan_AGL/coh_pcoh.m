cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
% coherence 12x2x12x5x448x448 (session x subject x trial x frequency x source x source)
load('coh_all.mat'); % values from 0-1
load('Pcoh_boolean.mat'); % partial coherence 12x2x12x5x448x448 (session x subject x trial x frequency x source x source
% session 1-12: in time sequence [20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;2022100402;20221005]
% odd number sessions are synch; even number sessions are synco
% subject 1: subject L; subject 2: subject R (the paring subject)
% trial 1-3: uncouple;   trial 4-6: L-lead;    trial 7-9: R-Lead;      trial 10-12: mutual;
% frequency 1-5: {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'}
% source 1-448: all cortical ROIs from Right to the left hemisphere

% select a trial and frequency to examine coherence and partial coherence
ses=12;
subj=1;
tr=2;
freq=2;
% plot
subplot(121)
coh=squeeze(coh_all(ses,subj,tr,freq,:,:));
imagesc(coh);colorbar;colormap('jet');
clim([0 1]);
title('coh')
subplot(122)
Pcoh=squeeze(Pcoh_boolean(ses,subj,tr,freq,:,:));
imagesc(Pcoh);colorbar;colormap('jet');
clim([0 1]);
title('Pcoh boolean')
sgtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr) ' freq ' num2str(freq)]);


%% Examine edges from Partial coherence in different conditions
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('SC.mat')
n_in=sum(triu(SC,1),'all'); % total number of edges inside SC
n_out=sum(triu(~SC,1),'all'); % total number of edges outside SC

% compute number of edges 
NedgeIn_coh=nan(12,2,12,5); % inside SC
NedgeOut_coh=nan(12,2,12,5); % outside SC
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                G=squeeze(Pcoh_boolean(ses,subj,tr,freq,:,:));
                NedgeIn_coh(ses,subj,tr,freq) = sum(G.*triu(SC,1),'all');
                NedgeOut_coh(ses,subj,tr,freq) =  sum(G.*triu(~SC,1),'all');
            end
        end
    end
end

% organzied into 2 syn type and 4 conditons
NedgeIn_coh4=cell(2,4,5); 
NedgeOut_coh4=cell(2,4,5);
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                if ismember(tr,[1 2 3]) && ismember(ses,[1:2:11]); % uncouple synch
                    NedgeIn_coh4{1,1,freq}=[NedgeIn_coh4{1,1,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{1,1,freq}=[NedgeOut_coh4{1,1,freq} NedgeOut_coh(ses,subj,tr,freq)];
                elseif ismember(tr,[1 2 3]) && ismember(ses,[2:2:12]); % uncouple synco
                    NedgeIn_coh4{2,1,freq}=[NedgeIn_coh4{2,1,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{2,1,freq}=[NedgeOut_coh4{2,1,freq} NedgeOut_coh(ses,subj,tr,freq)];

                elseif ((ismember(tr,[4:6]) && subj==1) || (ismember(tr,[7:9]) && subj==2)) && ismember(ses,[1:2:11]); % leading synch
                    NedgeIn_coh4{1,2,freq}=[NedgeIn_coh4{1,2,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{1,2,freq}=[NedgeOut_coh4{1,2,freq} NedgeOut_coh(ses,subj,tr,freq)];
                elseif ((ismember(tr,[4:6]) && subj==1) || (ismember(tr,[7:9]) && subj==2)) && ismember(ses,[2:2:12]); % leading synco
                    NedgeIn_coh4{2,2,freq}=[NedgeIn_coh4{2,2,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{2,2,freq}=[NedgeOut_coh4{2,2,freq} NedgeOut_coh(ses,subj,tr,freq)];

                elseif ((ismember(tr,[4:6]) && subj==2) || (ismember(tr,[7:9]) && subj==1)) && ismember(ses,[1:2:11]); % following synch
                    NedgeIn_coh4{1,3,freq}=[NedgeIn_coh4{1,3,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{1,3,freq}=[NedgeOut_coh4{1,3,freq} NedgeOut_coh(ses,subj,tr,freq)];
                elseif ((ismember(tr,[4:6]) && subj==2) || (ismember(tr,[7:9]) && subj==1)) && ismember(ses,[2:2:12]); % following synco
                    NedgeIn_coh4{2,3,freq}=[NedgeIn_coh4{2,3,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{2,3,freq}=[NedgeOut_coh4{2,3,freq} NedgeOut_coh(ses,subj,tr,freq)];

                elseif ismember(tr,[10:12]) && ismember(ses,[1:2:11]); % mutual synch
                    NedgeIn_coh4{1,4,freq}=[NedgeIn_coh4{1,4,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{1,4,freq}=[NedgeOut_coh4{1,4,freq} NedgeOut_coh(ses,subj,tr,freq)];
                else ismember(tr,[10:12]) && ismember(ses,[2:2:12]); % mutual synco
                    NedgeIn_coh4{2,4,freq}=[NedgeIn_coh4{2,4,freq} NedgeIn_coh(ses,subj,tr,freq)];
                    NedgeOut_coh4{2,4,freq}=[NedgeOut_coh4{2,4,freq} NedgeOut_coh(ses,subj,tr,freq)];
                end
            end
        end
    end
end

% compute the mean and standard error
NedgeIn_coh4mean=nan(2,4,5);
NedgeOut_coh4mean=nan(2,4,5);
NedgeIn_coh4ste=nan(2,4,5);
NedgeOut_coh4ste=nan(2,4,5);
for syn=1:2
    for con=1:4
        for freq=1:5
            NedgeIn_coh4mean(syn,con,freq)=mean(NedgeIn_coh4{syn,con,freq});
            NedgeOut_coh4mean(syn,con,freq)=mean(NedgeOut_coh4{syn,con,freq});
            NedgeIn_coh4ste(syn,con,freq)=std(NedgeIn_coh4{syn,con,freq})/sqrt(length(NedgeIn_coh4{syn,con,freq}));
            NedgeOut_coh4ste(syn,con,freq)=std(NedgeOut_coh4{syn,con,freq})/sqrt(length(NedgeOut_coh4{syn,con,freq}));
        end
    end
end

%% plotting
% load labels and colors for plots
run plotting_scheme.m
figure
for syn =1:2
    subplot(1,2,syn);
    model_series=(squeeze(NedgeIn_coh4mean(syn,:,:)))'/n_in;
    model_error=(squeeze(NedgeIn_coh4ste(syn,:,:)))'/n_in;
    b=bar(model_series,'grouped');
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(model_series);
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        b(i).FaceColor=condicolors(i,:);
        x(i,:) = b(i).XEndPoints;
    end
    % Plot the errorbars 
    hold on;
    errorbar(x', model_series, model_error,'k','linestyle','none','LineWidth',2);
    ylabel(['prt of edges'])
    xticks([1:5])
    xticklabels(bandlabels)
    
    legend(condi4names)
    ylim([0.05 0.3])
    title(['inside SC: ' syn2names{syn}],'color',syn2colors(syn,:));
end
