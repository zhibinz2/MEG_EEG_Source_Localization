clear
%% cov -> coh
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
% aggreate coh and pcoh inside SC for 2x4 condi
sc_coh=cell(2,4,5); % 32x32x36 trial within each cell
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                coh=squeeze(coh32_all_sorted(ses,subj,tr,freq,:,:)); 

                if ismember(tr,[1 2 3]) && ismember(ses,[1:2:11]); % uncouple synch
                    sc_coh{1,1,freq}=cat(3,sc_coh{1,1,freq},coh);
                    
                elseif ismember(tr,[1 2 3]) && ismember(ses,[2:2:12]); % uncouple synco
                    sc_coh{2,1,freq}=cat(3,sc_coh{2,1,freq},coh);


                elseif ((ismember(tr,[4:6]) && subj==1) || (ismember(tr,[7:9]) && subj==2)) && ismember(ses,[1:2:11]); % leading synch
                    sc_coh{1,2,freq}=cat(3,sc_coh{1,2,freq},coh);
                    
                elseif ((ismember(tr,[4:6]) && subj==1) || (ismember(tr,[7:9]) && subj==2)) && ismember(ses,[2:2:12]); % leading synco
                    sc_coh{2,2,freq}=cat(3,sc_coh{2,2,freq},coh);
                    

                elseif ((ismember(tr,[4:6]) && subj==2) || (ismember(tr,[7:9]) && subj==1)) && ismember(ses,[1:2:11]); % following synch
                    sc_coh{1,3,freq}=cat(3,sc_coh{1,3,freq},coh);
                    
                elseif ((ismember(tr,[4:6]) && subj==2) || (ismember(tr,[7:9]) && subj==1)) && ismember(ses,[2:2:12]); % following synco
                    sc_coh{2,3,freq}=cat(3,sc_coh{2,3,freq},coh);
                    

                elseif ismember(tr,[10:12]) && ismember(ses,[1:2:11]); % mutual synch
                    sc_coh{1,4,freq}=cat(3,sc_coh{1,4,freq},coh);
                    
                else ismember(tr,[10:12]) && ismember(ses,[2:2:12]); % mutual synco
                    sc_coh{2,4,freq}=cat(3,sc_coh{2,4,freq},coh);
                    
                end
            end
        end
    end
end
% paring states
% sc_coh=cell(2,4,5); % 32x32x36 trial within each cell
corr_coh_sampl=cell(3,5,2,32,32);% 3 states x 5 freq x 2 subjects x 32x32 x 36 trial within each cell
for freq =1:5
    for i=1:32
        for j=1:32
        % uncouple
            corr_coh_sampl{1,freq,1,i,j}=squeeze(sc_coh{1,1,freq}(i,j,:)); % subject L
            corr_coh_sampl{1,freq,2,i,j}=squeeze(sc_coh{2,1,freq}(i,j,:)); % subject R
    
            % unidirectional
            corr_coh_sampl{2,freq,1,i,j}=cat(1,squeeze(sc_coh{1,2,freq}(i,j,:)),squeeze(sc_coh{2,2,freq}(i,j,:))); % Leader
            corr_coh_sampl{2,freq,2,i,j}=cat(1,squeeze(sc_coh{1,3,freq}(i,j,:)),squeeze(sc_coh{2,3,freq}(i,j,:))); % Follower
    
            % bidirectional
            corr_coh_sampl{3,freq,1,i,j}=squeeze(sc_coh{1,4,freq}(i,j,:)); % subject L
            corr_coh_sampl{3,freq,2,i,j}=squeeze(sc_coh{2,4,freq}(i,j,:)); % subject R

        end
    end
end
% compute correlation
corr_coh=nan(3,5,32,32); % 3 states x 5 freq x 32 x 32 edges
for st=1:3
    for freq=1:5
        for i=1:32
            for j=1:32
                A=corr_coh_sampl{st,freq,1,i,j};
                B=corr_coh_sampl{st,freq,2,i,j};
                corr_coh(st,freq,i,j)=corr(A,B);
            end
        end
    end
end
% find the sources with correlations of 'Independent' < 'Unidirectional' < 'Bidirectional'
% corr_matched=nan(3,5,32,32);
corr_matched3=zeros(5,32,32);
for freq=1:5
    for i=1:32
        for j=1:32
%            if ((corr_coh(1,freq,i,j)+0.2) < corr_coh(2,freq,i,j)) && (corr_coh(2,freq,i,j) < (corr_coh(3,freq,i,j)-0.2)) ...
%                     && (abs(corr_coh(1,freq,i,j))<0.4) && (corr_coh(3,freq,i,j)>0.6)
           if (abs(corr_coh(1,freq,i,j))<0.4) && (corr_coh(3,freq,i,j)>0.6) && (corr_coh(2,freq,i,j)>0.4) ...
               && (corr_coh(2,freq,i,j) < (corr_coh(3,freq,i,j)))
%                corr_matched(1,freq,i,j)=corr_coh(1,freq,i,j);
%                corr_matched(2,freq,i,j)=corr_coh(2,freq,i,j);
%                corr_matched(3,freq,i,j)=corr_coh(3,freq,i,j);
                corr_matched3(freq,i,j)=1;
           end
        end
    end
end
figure
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/Coh_networkx
run plotting_scheme.m
direction3names={'Independent','Unidirectional','Bidirectional'};
dire3colors=[darkgreen;brown;megenta];
clf
for freq=1:5
    subplot(1,5,freq)
%     hold on
%     plot(1:32,diag(squeeze(corr_matched(1,freq,:,:))),'.','color',dire3colors(1,:),'MarkerSize',12);
%     plot(1:32,diag(squeeze(corr_matched(2,freq,:,:))),'.','color',dire3colors(2,:),'MarkerSize',12);
%     plot(1:32,diag(squeeze(corr_matched(3,freq,:,:))),'.','color',dire3colors(3,:),'MarkerSize',12);
%     xlabel('channels');ylabel('correlation');
%     ylim([-0.4 1]);
%     title(bandlabels{freq});
%     hold off;
%     if freq ==1
%         legend(direction3names)
%     end
%     grid on
    imagesc(squeeze(corr_matched3(freq,:,:)));colorbar;clim([0 1]);title(bandlabels{freq});
    xlabel('channel');ylabel('channel');
end
sgtitle({'Correlation of coherence between subject pair', ...
    'show only sources with correlation of Independent<Unidirectional<Bidirectional && abs(Independent)<0.4 && Unidirectional>0.4 && Bidirectional>0.6'});
set(gcf,'color','w'); 

%% pcov -> pcoh
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
