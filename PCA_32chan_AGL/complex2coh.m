%% color scheme
red   = [1 0 0];
pink  = [1 0.65 0.75];
black = [0 0 0];
white = [1 1 1];
blue  = [0 0 1];
mediumblue = [0 0.4 0.7];
green = [0 1 0];
darkgreen = [0 0.5 0];
grey  = [0.5 0.5 0.5];
yellow  = [1 1 0];
deepyellow  = [1 0.8 0.2];
gold = [212/255 175/255 55/255];
brown = [150/255 75/255 0];
megenta = [1 0 1];
cyan = [0 1 1]; 
purple = [0.6 0.1 0.9];
matlab_blue=[0 0.4470 0.7410];
matlab_orange=[0.8500 0.3250 0.0980];
matlab_gold=[0.9290 0.6940 0.1250];
matlab_purple=[0.4940 0.1840 0.5560];
matlab_green=[0.4660 0.6740 0.1880];
matlab_cyan=[0.3010 0.7450 0.9330];
matlab_red=[0.6350 0.0780 0.1840];
% combine colors
freq5colors=[matlab_blue;matlab_cyan;matlab_green;matlab_orange;matlab_red];
condicolors=[darkgreen;red;blue;megenta;purple;purple];
dire3colors=[darkgreen;brown;megenta];
syn2colors=[darkgreen;pink];
HNLcolors = [darkgreen; deepyellow; pink];



%%
open /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util/normalizeCSD.m
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('Complex_X_op3_all.mat');

ses=1;
subj=1;
tr=11;
freq=5;
cov_c=squeeze(Complex_X_op3_all(ses,subj,tr,freq,:,:));
% Complex_X_op3_all.mat is 448x448 complex covariance after AGL
% condition sorted already
coh=normalizeCSD(cov_c);
imagesc(coh);colorbar;
vlim=0.005;
clim([0 vlim]);
title('coh')
subtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr) ' freq ' num2str(freq)]);

% normalize Complex_X_op3_all to real value partial coherence
% examine the number edges (using output from r2c in hilbert2cov.m)
NedgeIn_coh=nan(12,2,12,5);
NedgeOut_coh=nan(12,2,12,5);
Pcoh_all=nan(12,2,12,5,448,448);
tic
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                cov_c=squeeze(Complex_X_op3_all(ses,subj,tr,freq,:,:));
                newG=normalizeCSD(cov_c);
                Pcoh_all(ses,subj,tr,freq,:,:)=newG;
                NedgeIn_coh(ses,subj,tr,freq) = sum(sum(logical(newG).*triu(SC,1)));
                NedgeOut_coh(ses,subj,tr,freq) =  sum(sum(logical(newG).*triu(~SC,1)));
            end
        end
    end
end
toc % 12S
% cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
% save('Pcoh_all.mat','Pcoh_all','-v7.3');
% Pcoh_all.mat is real value partial coherence after AGL (448X448)
% with condition sorted already

ses=1;
subj=1;
tr=11;
freq=5;
Pcoh=squeeze(Pcoh_all(ses,subj,tr,freq,:,:));
imagesc(Pcoh);colorbar;
vlim=0.005;
clim([0 vlim]);
title('P_coh')
subtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr) ' freq ' num2str(freq)]);

% convert Pcoh_all to boolean
Pcoh_boolean=logical(Pcoh_all);
save('Pcoh_boolean.mat','Pcoh_boolean','-v7.3');
% Pcoh_boolean.mat is logical value partial coherence after AGL (448X448)
% with condition sorted already

ses=3;
subj=2;
tr=2;
freq=4;
Pcoh=squeeze(Pcoh_boolean(ses,subj,tr,freq,:,:));
imagesc(Pcoh);colorbar;
vlim=0.005;
clim([0 vlim]);
title('P_coh')
subtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr) ' freq ' num2str(freq)]);

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
                else ismember(tr,[10:12]); % mutual
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
NedgeIn_coh4ste=nan(4,5);
NedgeOut_coh4ste=nan(4,5);
for con=1:4
    for freq=1:5
        NedgeIn_coh4mean(con,freq)=mean(NedgeIn_coh4{con,freq});
        NedgeOut_coh4mean(con,freq)=mean(NedgeOut_coh4{con,freq});
        NedgeIn_coh4ste(con,freq)=std(NedgeIn_coh4{con,freq})/sqrt(length(NedgeIn_coh4{con,freq}));
        NedgeOut_coh4ste(con,freq)=std(NedgeOut_coh4{con,freq})/sqrt(length(NedgeOut_coh4{con,freq}));
    end
end

n_in=sum(triu(SC,1),'all');
n_out=sum(triu(~SC,1),'all');

figure;
for con=1:4
    subplot(1,4,con)
    hold on;
    plot(1:5,NedgeIn_coh4mean(con,:)/n_in,'g');
    plot(1:5,NedgeOut_coh4mean(con,:)/n_out,'r');
    ylabel('n edges (prt)')
    xlabel('frequency band')
    title(['condition ' num2str(con)])
    subtitle(['percentage of edges in and outside SC']);
    xticks([1:5])
    xticklabels(bandlabels)
    ylim([0 0.2])
    hold off;
end

figure;
b=bar(NedgeIn_coh4mean'/n_in,'grouped');
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(NedgeIn_coh4mean'/n_in);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    b(i).FaceColor=condicolors(i,:);
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars 
hold on;
errorbar(x',NedgeIn_coh4mean'/n_in,NedgeIn_coh4ste'/n_in,'k','linestyle','none','LineWidth',2);
title('inside SC')
ylabel(['prt of edges'])
xticks([1:5])
xticklabels(bandlabels)
condi4names={'Uncouple','Leading','Following','Mutual'};
legend(condi4names)
ylim([0.1 0.25])

%% sort through condtions and syn-types (only inside SC)
% boolean matrix average
NedgeIn_coh4=cell(2,4,5);
NedgeOut_coh4=cell(2,4,5);
tic
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
toc

bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};

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

n_in=sum(triu(SC,1),'all');
n_out=sum(triu(~SC,1),'all');

syn2names={'syncH','syncO'};
syn=1;
syn=2;

figure;
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
condi4names={'Uncouple','Leading','Following','Mutual'};
legend(condi4names)
ylim([0.05 0.3])
title(['inside SC: ' syn2names{syn}]);

%% organize for networkx both coh and pcoh

%% correlation with granger





