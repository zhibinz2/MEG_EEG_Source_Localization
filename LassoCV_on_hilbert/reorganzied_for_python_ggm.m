cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('hilbert_dataCov_all.mat');

%% preapare
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/
for ses=12%12
    tic
    for subj=1:2
        for tr=1:12
            % condition not sorted
            load(['./cortical_source_data/' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat']);

            hilbert_dataCov=cell(1,5);
            for freq=1:5
                
                % resample
                downsample_data=resample(double(agr_source_data),1,downsample,'Dimension',1);
                filterd_data = filter(filt_ds{freq},downsample_data);
                hilbertdata = hilbert(filterd_data'); 

                sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));
                sourceDataReal = sourceDataReal*(1/mean(abs(sourceDataReal(:)))); % normalize data

                hilbert_dataCov{freq} = cov(sourceDataReal');

                save(['./hilbert_datacov/' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat'],'hilbert_dataCov');
            end
        end
    end
    toc
end

%% Combine the precision matrix after applying lasso in python on multiple computers
clear
tic
cd /ssd/zhibin/Cleaned_sourcedata/cortical_source_data/python_lasso
prec_all_ses_1_12=zeros(12,2,12,5,896,896);

clear prec_all
load('prec_all_ses_0_1.mat')
prec_all_ses_1_2=prec_all;
ses=2;subj=2;tr=5;freq=1;
prec=squeeze(prec_all_ses_1_2(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(1:2,:,:,:,:,:)=prec_all_ses_1_2(1:2,:,:,:,:,:);
clear prec_all_ses_1_2 prec_all

load('prec_all_ses_2_3.mat') % hnlb
prec_all_ses_3_4=prec_all;
ses=4;subj=2;tr=6;freq=5;
prec=squeeze(prec_all_ses_3_4(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(3:4,:,:,:,:,:)=prec_all_ses_3_4(3:4,:,:,:,:,:);
clear prec_all_ses_3_4 prec_all

load('prec_all_ses_4.mat') % ramesh
prec_all_ses_5=prec_all;
ses=10;subj=1;tr=6;freq=2;
prec=squeeze(prec_all_ses_5(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(5,:,:,:,:,:)=prec_all_ses_5(5,:,:,:,:,:);
clear prec_all_ses_5 prec_all

load('prec_all_ses_5.mat') % hnla
prec_all_ses_6=prec_all;
ses=6;subj=2;tr=6;freq=1;
prec=squeeze(prec_all_ses_6(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(6,:,:,:,:,:)=prec_all_ses_6(6,:,:,:,:,:);
clear prec_all_ses_6 prec_all

load('prec_all_ses_6_7.mat');
prec_all_ses_7_8=prec_all;
ses=7;subj=1;tr=5;freq=3;
prec=squeeze(prec_all_ses_7_8(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(7:8,:,:,:,:,:)=prec_all_ses_7_8(7:8,:,:,:,:,:);
clear prec_all_ses_7_8 prec_all

load('prec_all_ses_8.mat') 
prec_all_ses_9=prec_all;
ses=9;subj=1;tr=6;freq=2;
prec=squeeze(prec_all_ses_9(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(9,:,:,:,:,:)=prec_all_ses_9(9,:,:,:,:,:);
clear prec_all_ses_9 prec_all

load('prec_all_ses_9.mat') % hnla
prec_all_ses_10=prec_all;
ses=10;subj=1;tr=6;freq=2;
prec=squeeze(prec_all_ses_10(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(10,:,:,:,:,:)=prec_all_ses_10(10,:,:,:,:,:);
clear prec_all_ses_10 prec_all

load('prec_all_ses_10.mat') % hnla works % ramesh?
prec_all_ses_11=prec_all;
ses=11;subj=2;tr=3;freq=3;
prec=squeeze(prec_all_ses_11(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(11,:,:,:,:,:)=prec_all_ses_11(11,:,:,:,:,:);
clear prec_all_ses_11 prec_all

load('prec_all_ses_11.mat') % ramesh
prec_all_ses_12=prec_all;
ses=12;subj=2;tr=5;freq=2;
prec=squeeze(prec_all_ses_12(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(12,:,:,:,:,:)=prec_all_ses_12(12,:,:,:,:,:);
clear prec_all_ses_12 prec_all

cd /ssd/zhibin/1overf/Cleaned_sourcedata/cortical_source_data/python_lasso
save('prec_all_ses_1_12.mat','prec_all_ses_1_12','-v7.3'); % condition not sorted
toc % 300s

figure
sumpre_all=nan(12,2,12,5);
for ses=1:12 
    for subj=1:2
        for tr=1:12
            for freq=1:5
                prec=squeeze(prec_all_ses_1_12(ses,subj,tr,freq,:,:));
                imagesc(logical(prec));colorbar;
                sumprec=sum(triu(prec,1),"all");
                title(['sum of prec = ' num2str(sumprec)]);
                sumpre_all(ses,subj,tr,freq)=sumprec;
                subtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' tr ' num2str(tr) ...
                    ' freq ' num2str(freq)])
                pause(0.1)
            end
        end
    end
end

%%
% convert to complex prec
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
Complex_prec_all=nan(12,2,12,5,448,448); % condition not sorted
tic
for ses=1:12
    for subj=1:2
        for tr=1:12
            for freq=1:5
                Complex_prec_all(ses,subj,tr,freq,:,:)=r2c(squeeze(prec_all_ses_1_12(ses,subj,tr,freq,:,:)));
            end
        end
    end
end
toc

% normalize to real coherence
cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
coh_lasso=nan(12,2,12,5,448,448);
tic
for ses=1:12
    for subj=1:2
        for tr=1:12
            for freq=1:5
                coh_lasso(ses,subj,tr,freq,:,:)=normalizeCSD(squeeze(Complex_prec_all(ses,subj,tr,freq,:,:)));
            end
        end
    end
end
toc % 10s
cd /ssd/zhibin/1overf/Cleaned_sourcedata/cortical_source_data/python_lasso
Pcoh_lasso=coh_lasso;
save('Pcoh_lasso.mat','Pcoh_lasso','-v7.3');
% condition not sorted


figure
for ses=1:12 
    for subj=1:2
        for tr=1:12
            for freq=1:5
                matshow=logical(squeeze(coh_lasso(ses,subj,tr,freq,:,:)));
                imagesc(matshow);
                sumedge=sum(triu(matshow,1),"all");
                title(['sum of prec = ' num2str(sumedge)]);
                sumedge_all(ses,subj,tr,freq)=sumedge;
                subtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' tr ' num2str(tr) ...
                    ' freq ' num2str(freq)])
                pause(0.1)
            end
        end
    end
end

%% compute # edges in and out
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
                G=logical(squeeze(coh_lasso(ses,subj,tr,freq,:,:)));
                NedgeIn_coh(ses,subj,tr,freq) = sum(G.*triu(SC,1),'all');
                NedgeOut_coh(ses,subj,tr,freq) =  sum(G.*triu(~SC,1),'all');
            end
        end
    end
end

%% sort condition here


%% 2X4 condi (wrong, condition not sorted)
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
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/Coh_networkx
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
    % ylim([0.05 0.3])
    title(['inside SC: ' syn2names{syn}],'color',syn2colors(syn,:));
end


%% Examine average coh and pcoh in 2x4 conditions inside SC 
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('SC.mat');
% aggreate coh and pcoh inside SC for 2x4 condi
sc_coh_lasso=cell(2,4,5); % 448x448x36 trial within each cell
tic
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                coh_las=squeeze(coh_lasso(ses,subj,tr,freq,:,:)).*SC; % investigate only within SC

                if ismember(tr,[1 2 3]) && ismember(ses,[1:2:11]); % uncouple synch
                    sc_coh_lasso{1,1,freq}=cat(3,sc_coh_lasso{1,1,freq},coh_las);
                elseif ismember(tr,[1 2 3]) && ismember(ses,[2:2:12]); % uncouple synco
                    sc_coh_lasso{2,1,freq}=cat(3,sc_coh_lasso{2,1,freq},coh_las);

                elseif ((ismember(tr,[4:6]) && subj==1) || (ismember(tr,[7:9]) && subj==2)) && ismember(ses,[1:2:11]); % leading synch
                    sc_coh_lasso{1,2,freq}=cat(3,sc_coh_lasso{1,2,freq},coh_las);
                elseif ((ismember(tr,[4:6]) && subj==1) || (ismember(tr,[7:9]) && subj==2)) && ismember(ses,[2:2:12]); % leading synco
                    sc_coh_lasso{2,2,freq}=cat(3,sc_coh_lasso{2,2,freq},coh_las);

                elseif ((ismember(tr,[4:6]) && subj==2) || (ismember(tr,[7:9]) && subj==1)) && ismember(ses,[1:2:11]); % following synch
                    sc_coh_lasso{1,3,freq}=cat(3,sc_coh_lasso{1,3,freq},coh_las);
                elseif ((ismember(tr,[4:6]) && subj==2) || (ismember(tr,[7:9]) && subj==1)) && ismember(ses,[2:2:12]); % following synco
                    sc_coh_lasso{2,3,freq}=cat(3,sc_coh_lasso{2,3,freq},coh_las);

                elseif ismember(tr,[10:12]) && ismember(ses,[1:2:11]); % mutual synch
                    sc_coh_lasso{1,4,freq}=cat(3,sc_coh_lasso{1,4,freq},coh_las);
                else ismember(tr,[10:12]) && ismember(ses,[2:2:12]); % mutual synco
                    sc_coh_lasso{2,4,freq}=cat(3,sc_coh_lasso{2,4,freq},coh_las);
                end
            end
        end
    end
end
toc % 22s

% mean coh
m_pcoh=zeros(2,4,5,448,448);
for syn = 1:2
    for condi = 1:4
        for freq = 1:5
            m_pcoh(syn,condi,freq,:,:)=mean(sc_coh_lasso{syn,condi,freq},3);
        end
    end
end

syn=2;
condi=4;
freq=5;
imagesc(logical(squeeze(m_pcoh(syn,condi,freq,:,:))));

% examine the mean in a plot
A=m_pcoh;
% plot
for freq=1:5
    figure('units','normalized','outerposition',[0 0 1 1]);
    sgtitle([bandlabels{freq} ' pcoh'])
    for syn =1:2
        for condi=1:4
            subplot(2,4,4*(syn-1)+condi)
            imagesc(logical(squeeze(A(syn,condi,freq,:,:))));colorbar;
            colormap('hot');
            title([syn2names{syn} ' ' condi4names{condi}]);
        end
    end
end
% plot
for freq=1:5
    figure('units','normalized','outerposition',[0 0 1 1]);
    sgtitle([bandlabels{freq} ' pcoh'])
    for syn =1:2
        for condi=1:4
            subplot(2,4,4*(syn-1)+condi)
            imagesc(squeeze(A(syn,condi,freq,:,:)));colorbar;
            colormap('jet');
            clim([-0.001 0.001])
            title([syn2names{syn} ' ' condi4names{condi}]);
        end
    end
end


%% corrrelation paring condition
% sc_coh=cell(2,4,5); % 448x448x36 trial within each cell
% corr_coh_lasso_sampl=cell(2,4,5); % 448x448x36 trial within each cell
corr_coh_lasso_sampl=cell(3,5,2,448,448);% 3 states x 5 freq x 2 subjects x 448x448 x 36 trial within each cell

for freq =1:5
    tic
    for i=1:448
        for j=1:448
        % uncouple
        corr_pcoh_sampl{1,freq,1,i,j}=squeeze(coh_lasso{1,1,freq}(i,j,:)); % subject L
        corr_pcoh_sampl{1,freq,2,i,j}=squeeze(coh_lasso{2,1,freq}(i,j,:)); % subject R
        % unidirectional
        corr_coh_sampl{2,freq,1,i,j}=cat(1,squeeze(sc_coh{1,2,freq}(i,j,:)),squeeze(sc_coh{2,2,freq}(i,j,:))); % Leader
        corr_coh_sampl{2,freq,2,i,j}=cat(1,squeeze(sc_coh{1,3,freq}(i,j,:)),squeeze(sc_coh{2,3,freq}(i,j,:))); % Follower
        corr_pcoh_sampl{2,freq,1,i,j}=cat(1,squeeze(sc_coh_lasso{1,2,freq}(i,j,:)),squeeze(sc_coh_lasso{2,2,freq}(i,j,:))); % Leader
        corr_pcoh_sampl{2,freq,2,i,j}=cat(1,squeeze(sc_coh_lasso{1,3,freq}(i,j,:)),squeeze(sc_coh_lasso{2,3,freq}(i,j,:))); % Follower
        % bidirectional
        corr_coh_sampl{3,freq,1,i,j}=squeeze(sc_coh{1,4,freq}(i,j,:)); % subject L
        corr_coh_sampl{3,freq,2,i,j}=squeeze(sc_coh{2,4,freq}(i,j,:)); % subject R
        corr_pcoh_sampl{3,freq,1,i,j}=squeeze(sc_coh_lasso{1,4,freq}(i,j,:)); % subject L
        corr_pcoh_sampl{3,freq,2,i,j}=squeeze(sc_coh_lasso{2,4,freq}(i,j,:)); % subject R
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