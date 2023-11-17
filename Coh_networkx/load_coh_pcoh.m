cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
% coherence 12x2x12x5x448x448 (session x subject x trial x frequency x source x source)
load('coh_all.mat'); % values from 0-1 (12x2x12x5x448x448 : session x subject x trial x frequency x source x source)
load('Pcoh_boolean.mat'); % partial coherence 12x2x12x5x448x448 (session x subject x trial x frequency x source x source
% session 1-12: in time sequence [20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;2022100402;20221005]
% odd number sessions are synch; even number sessions are synco
% subject 1: subject L; subject 2: subject R (the paring subject)
% trial 1-3: uncouple;   trial 4-6: L-lead;    trial 7-9: R-Lead;      trial 10-12: mutual;
% frequency 1-5: {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'}
% source 1-448: all cortical ROIs from Right to the left hemisphere

% select a trial and frequency to examine coherence and partial coherence
ses=2;
subj=1;
tr=1;
freq=2;
% plot
subplot(121)
coh=squeeze(coh_all(ses,subj,tr,freq,:,:));
% imagesc(coh);colorbar;colormap('jet');
% imagesc(SC);colorbar;colormap('jet');
imagesc(coh.*SC);colorbar;colormap('jet');
clim([0 1]);
title('coh')
subplot(122)
Pcoh=squeeze(Pcoh_boolean(ses,subj,tr,freq,:,:));
imagesc(Pcoh.*SC);colorbar;colormap('jet');
clim([0 1]);
title('Pcoh boolean')
sgtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr) ' freq ' num2str(freq)]);

%% Examine average coh and pcoh in 2x4 conditions inside SC 
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('SC.mat');
% aggreate coh and pcoh inside SC for 2x4 condi
sc_coh=cell(2,4,5); % 448x448x36 trial within each cell
sc_pcoh=cell(2,4,5); % 448x448x36 trial within each cell
tic
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                coh=squeeze(coh_all(ses,subj,tr,freq,:,:)).*SC; % investigate only within SC
                pcoh=squeeze(Pcoh_boolean(ses,subj,tr,freq,:,:)).*SC;

                if ismember(tr,[1 2 3]) && ismember(ses,[1:2:11]); % uncouple synch
                    sc_coh{1,1,freq}=cat(3,sc_coh{1,1,freq},coh);
                    sc_pcoh{1,1,freq}=cat(3,sc_pcoh{1,1,freq},pcoh);
                elseif ismember(tr,[1 2 3]) && ismember(ses,[2:2:12]); % uncouple synco
                    sc_coh{2,1,freq}=cat(3,sc_coh{2,1,freq},coh);
                    sc_pcoh{2,1,freq}=cat(3,sc_pcoh{2,1,freq},pcoh);

                elseif ((ismember(tr,[4:6]) && subj==1) || (ismember(tr,[7:9]) && subj==2)) && ismember(ses,[1:2:11]); % leading synch
                    sc_coh{1,2,freq}=cat(3,sc_coh{1,2,freq},coh);
                    sc_pcoh{1,2,freq}=cat(3,sc_pcoh{1,2,freq},pcoh);
                elseif ((ismember(tr,[4:6]) && subj==1) || (ismember(tr,[7:9]) && subj==2)) && ismember(ses,[2:2:12]); % leading synco
                    sc_coh{2,2,freq}=cat(3,sc_coh{2,2,freq},coh);
                    sc_pcoh{2,2,freq}=cat(3,sc_pcoh{2,2,freq},pcoh);

                elseif ((ismember(tr,[4:6]) && subj==2) || (ismember(tr,[7:9]) && subj==1)) && ismember(ses,[1:2:11]); % following synch
                    sc_coh{1,3,freq}=cat(3,sc_coh{1,3,freq},coh);
                    sc_pcoh{1,3,freq}=cat(3,sc_pcoh{1,3,freq},pcoh);
                elseif ((ismember(tr,[4:6]) && subj==2) || (ismember(tr,[7:9]) && subj==1)) && ismember(ses,[2:2:12]); % following synco
                    sc_coh{2,3,freq}=cat(3,sc_coh{2,3,freq},coh);
                    sc_pcoh{2,3,freq}=cat(3,sc_pcoh{2,3,freq},pcoh);

                elseif ismember(tr,[10:12]) && ismember(ses,[1:2:11]); % mutual synch
                    sc_coh{1,4,freq}=cat(3,sc_coh{1,4,freq},coh);
                    sc_pcoh{1,4,freq}=cat(3,sc_pcoh{1,4,freq},pcoh);
                else ismember(tr,[10:12]) && ismember(ses,[2:2:12]); % mutual synco
                    sc_coh{2,4,freq}=cat(3,sc_coh{2,4,freq},coh);
                    sc_pcoh{2,4,freq}=cat(3,sc_pcoh{2,4,freq},pcoh);
                end
            end
        end
    end
end
toc % 44s

% mean coh and pcoh
m_coh=zeros(2,4,5,448,448);
m_pcoh=zeros(2,4,5,448,448);
for syn = 1:2
    for condi = 1:4
        for freq = 1:5
            m_coh(syn,condi,freq,:,:)=mean(sc_coh{syn,condi,freq},3);
            m_pcoh(syn,condi,freq,:,:)=mean(sc_pcoh{syn,condi,freq},3);
        end
    end
end

% examine the mean in a plot
% select coh or pcoh to plot
A=m_coh;
A=m_pcoh;
% plot
for freq=1:5
    figure;
    sgtitle([bandlabels{freq} ' pcoh'])
    for syn =1:2
        for condi=1:4
            subplot(2,4,4*(syn-1)+condi)
            imagesc(squeeze(A(syn,condi,freq,:,:)));colorbar;
            colormap('jet');
            title([syn2names{syn} ' ' condi4names{condi}]);
        end
    end
end

% The above look all very similar
% need zscoreing the 2x4 conditions
zr_coh=nan(8,5,448,448);% reorganzied into 8 x 5 x 488 x 488 then zscore along the first dimention
zr_pcoh=nan(8,5,448,448);
for freq =1:5
    i=1;
    for syn=1:2
        for condi=1:4
            zr_coh(i,freq,:,:)=squeeze(m_coh(syn,condi,freq,:,:));
            zr_pcoh(i,freq,:,:)=squeeze(m_pcoh(syn,condi,freq,:,:));
            i=i+1;
        end
    end    
end
% convert to z-score
z_coh=nan(8,5,448,448);% reorganzied into 8 x 5 x 488 x 488 then zscore along the first dimention
z_pcoh=nan(8,5,448,448);
for freq =1:5
    z_coh(:,freq,:,:)=zscore(zr_coh(:,freq,:,:),[],1);
    z_pcoh(:,freq,:,:)=zscore(zr_pcoh(:,freq,:,:),[],1);
end
% plot
A=z_coh;
A=z_pcoh;
for freq=1:5
    figure('units','normalized','outerposition',[0 0 1 1]);
    sgtitle([bandlabels{freq} '   zscore among the 8 condition.  pcoh'])
    i = 1;
    for syn =1:2
        for condi=1:4
            subplot(2,4,i)
            imagesc(squeeze(A(i,freq,:,:)));colorbar;
            colormap('hotncold');
            title([syn2names{syn} ' ' condi4names{condi}]);
            clim([-2.5 2.5])
            i=i+1;
        end
    end
end

% try difference to the mean
d_coh=nan(8,5,448,448);% reorganzied into 8 x 5 x 488 x 488 then zscore along the first dimention
d_pcoh=nan(8,5,448,448);
for freq =1:5
    d_coh(:,freq,:,:)=zr_coh(:,freq,:,:)-repmat(mean(zr_coh(:,freq,:,:),1),8,1,1,1);
    d_pcoh(:,freq,:,:)=zr_pcoh(:,freq,:,:)-repmat(mean(zr_pcoh(:,freq,:,:),1),8,1,1);
end
% plot
A=d_coh;
A=d_pcoh;
for freq=1:5
    figure('units','normalized','outerposition',[0 0 1 1]);
    sgtitle(bandlabels{freq})
    i = 1;
    for syn =1:2
        for condi=1:4
            subplot(2,4,i)
            imagesc(squeeze(A(i,freq,:,:)));colorbar;
            colormap('jet');
            title([syn2names{syn} ' ' condi4names{condi}]);
            i=i+1;
        end
    end
end



%% Examine number of edges from Partial coherence in different conditions
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('SC.mat');
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


%% corrrelation paring condition
% sc_coh=cell(2,4,5); % 448x448x36 trial within each cell
% sc_pcoh=cell(2,4,5); % 448x448x36 trial within each cell
corr_coh_sampl=cell(3,5,2,448,448);% 3 states x 5 freq x 2 subjects x 448x448 x 36 trial within each cell
corr_pcoh_sampl=cell(3,5,2,448,448);
for freq =1:5
    tic
    for i=1:448
        for j=1:448
        % uncouple
        corr_coh_sampl{1,freq,1,i,j}=squeeze(sc_coh{1,1,freq}(i,j,:)); % subject L
        corr_coh_sampl{1,freq,2,i,j}=squeeze(sc_coh{2,1,freq}(i,j,:)); % subject R
        corr_pcoh_sampl{1,freq,1,i,j}=squeeze(sc_pcoh{1,1,freq}(i,j,:)); % subject L
        corr_pcoh_sampl{1,freq,2,i,j}=squeeze(sc_pcoh{2,1,freq}(i,j,:)); % subject R
        % unidirectional
        corr_coh_sampl{2,freq,1,i,j}=cat(1,squeeze(sc_coh{1,2,freq}(i,j,:)),squeeze(sc_coh{2,2,freq}(i,j,:))); % Leader
        corr_coh_sampl{2,freq,2,i,j}=cat(1,squeeze(sc_coh{1,3,freq}(i,j,:)),squeeze(sc_coh{2,3,freq}(i,j,:))); % Follower
        corr_pcoh_sampl{2,freq,1,i,j}=cat(1,squeeze(sc_pcoh{1,2,freq}(i,j,:)),squeeze(sc_pcoh{2,2,freq}(i,j,:))); % Leader
        corr_pcoh_sampl{2,freq,2,i,j}=cat(1,squeeze(sc_pcoh{1,3,freq}(i,j,:)),squeeze(sc_pcoh{2,3,freq}(i,j,:))); % Follower
        % bidirectional
        corr_coh_sampl{3,freq,1,i,j}=squeeze(sc_coh{1,4,freq}(i,j,:)); % subject L
        corr_coh_sampl{3,freq,2,i,j}=squeeze(sc_coh{2,4,freq}(i,j,:)); % subject R
        corr_pcoh_sampl{3,freq,1,i,j}=squeeze(sc_pcoh{1,4,freq}(i,j,:)); % subject L
        corr_pcoh_sampl{3,freq,2,i,j}=squeeze(sc_pcoh{2,4,freq}(i,j,:)); % subject R
        end
    end
    toc % 15s 
end
% 1.25 min


corr_coh=nan(3,5,448,448); % 3 states x 5 freq x 448 x 448 edges
corr_pcoh=nan(3,5,448,448); % 3 states x 5 freq x 448 x 448 edges
for st=1:3
    for freq=1:5
        tic
        for i=1:448
            for j=1:448
                A=corr_coh_sampl{st,freq,1,i,j};
                B=corr_coh_sampl{st,freq,2,i,j};
                corr_coh(st,freq,i,j)=corr(A,B);
               
                A=corr_pcoh_sampl{st,freq,1,i,j};
                B=corr_pcoh_sampl{st,freq,2,i,j};
                corr_pcoh(st,freq,i,j)=corr(A,B);
            end
        end
        toc % 15 s
    end
end

% find the sources with correlations of 'Independent' < 'Unidirectional' < 'Bidirectional'
corr_dg_ctr_matched=nan(3,5,448);
for freq=1:5
    for sr=1:448
       if ((corr_dg_ctr(1,freq,sr)+0.2) < corr_dg_ctr(2,freq,sr)) && (corr_dg_ctr(2,freq,sr) < (corr_dg_ctr(3,freq,sr)-0.2)) ...
               && (abs(corr_dg_ctr(1,freq,sr))<0.4) && (corr_dg_ctr(3,freq,sr)>0.6)
           corr_dg_ctr_matched(1,freq,sr)=corr_dg_ctr(1,freq,sr);
           corr_dg_ctr_matched(2,freq,sr)=corr_dg_ctr(2,freq,sr);
           corr_dg_ctr_matched(3,freq,sr)=corr_dg_ctr(3,freq,sr);
       end
    end
end

figure
clf
direction3names={'Independent','Unidirectional','Bidirectional'};
dire3colors=[darkgreen;brown;megenta];
for freq=1:5
    subplot(5,1,freq)
    hold on
    plot(1:448,squeeze(corr_dg_ctr_matched(1,freq,:)),'.','color',dire3colors(1,:),'MarkerSize',12);
    plot(1:448,squeeze(corr_dg_ctr_matched(2,freq,:)),'.','color',dire3colors(2,:),'MarkerSize',12);
    plot(1:448,squeeze(corr_dg_ctr_matched(3,freq,:)),'.','color',dire3colors(3,:),'MarkerSize',12);
    xlabel('sources');ylabel('correlation');
    ylim([-0.4 1]);
    title(bandlabels{freq});
    hold off;
    if freq ==1
        legend(direction3names)
    end
    grid on
end
sgtitle({'Correlation of degreee centrality between subject pair', ...
    'show only sources with correlation of Independent+0.2<Unidirectional<Bidirectional-0.2 && abs(Independent)<0.4 && Bidirectional>0.6'});
set(gcf,'color','w'); 

%% compute #edges in and out
