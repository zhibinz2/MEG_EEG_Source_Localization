%% loop through all cortical source data to compute hilbert and cov in 5 freq
clear

cd /home/zhibinz2/Documents/GitHub/Cleaned_data/cortical_source_data
load('corti_ave_source_labl.mat')
load('corti_ave_source_coor.mat')
source_labels=corti_ave_source_labl{1,1,1};
source_coor=corti_ave_source_coor{1,1,1};

sampl_rate=2000;
srnew = 200;
downsample = 10;
passbands = [3.5 6.5; 7 10; 10.5 13.5; 14 20; 21 29];
bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
attenuation=60;
filt_ds=cell(1,length(passbands));
for freqBand=1:5
    % Select frequency
    passFreq1 = passbands(freqBand,1);
    passFreq2 = passbands(freqBand,2);
    filt_d = designfilt('bandpassiir','FilterOrder',20, ...
    'PassbandFrequency1',passFreq1,'PassbandFrequency2',passFreq2, ...
    'StopbandAttenuation1',attenuation,'PassbandRipple',0.2, ...
    'StopbandAttenuation2',attenuation,'SampleRate',srnew);
    filt_ds{freqBand}=filt_d;
end

seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
% for ses=1:12
%     mkdir(num2str(seeds(ses,:)));
% end
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/
for ses=12%12
    tic
    for subj=1:2
        for tr=1:12
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
% 1h 


%% load variables for penalty selection
clear
allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]); 
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
n_Lambdas=length(allLambdas); % number of lambda values
min_LamdaIn = min(allLambdas);

cd /home/zhibinz2/Documents/GitHub/Cleaned_data/cortical_source_data
load('corti_ave_source_labl.mat')
load('corti_ave_source_coor.mat')
source_labels=corti_ave_source_labl{1,1,1}; % every trial is the same, load one trial
source_coor=corti_ave_source_coor{1,1,1};

load('../../Virtual-Tractography/ForZhibin/processed_data/scale250_Connectome.mat');
SC=logical(fc(source_labels,source_labels));
sum(triu(SC,1),"all")


%% Organize hilber_dataCov_all
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
hilbert_dataCov_all=nan(12,2,12,5,896,896);
for ses=1:numSes
    tic
    clear conditions sortorders
    runid = num2str(seeds(ses,:));
    load(['../../Cleaned_data/clean_' runid '.mat'],'conditions');
    % sort order
    [x,sortorder]=sort(conditions);
    for subj=1:2
        for tr=1:12
            load(['../../Cleaned_data/hilbert_datacov/' ...
                    num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(sortorder(tr)) '.mat'] ...
                    ,'hilbert_dataCov');
            for freq=1:5
                hilbert_dataCov_all(ses,subj,tr,freq,:,:)=hilbert_dataCov{freq};
            end
        end
    end
    toc
end
% 1 min
% cd ../../Cleaned_data/hilbert_datacov/
% save('hilbert_dataCov_all.mat','hilbert_dataCov_all','-v7.3');


%% option 3 By subject
% We should also save the deviance values for the final fit, because they tell us how good a model it is.
% Option 3: Choose the penalties for each subject in each frequency band for all conditions.
% Each subject produces 24 trials, 3 in each condition.  So take 1 trial of each condition (8 trials 4 from each day) 
% and average covariance.  now you have 3 covariance matrices.  run the penalty selection function. .
% Go back and fit  each trial for that subject using that penalty.
% This way you would have to 12 x 5 penality optimizations.  But, I think this may be better, because the subject to subject differences may be greater than the conditon differences.
% Can you try this for the same frequency you are testing for option 1 in 1 subject.

% hilbert_dataCov_all=nan(12,2,12,5,894,894);
cd ../../Cleaned_data/hilbert_datacov
ave_hilcov_option3=nan(3,2,6,5,896,896); % 3 ensamble x 2 subjects x 6 double-sessions x 5 frequency x 894 x 894
for ensam = 1:3
    for subj =1:2
        for dl_ses=1:6
            for freq=1:5
                tic
                ave_hilcov_option3(ensam,subj,dl_ses,freq,:,:)=...
                    squeeze(mean(squeeze(mean(squeeze(hilbert_dataCov_all([1+2*(dl_ses-1):2+2*(dl_ses-1)], ...
                    subj,[1:3:12]+(ensam-1),freq,:,:)),2)),1));
                toc
            end
        end
    end
end

% penaltyselection for 1 subject
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
penalizationIn_op3=nan(2,6,5);
penalizationOut_op3=nan(2,6,5);
minDev_op3=nan(2,6,5);
parfor subj=1%:2
    for dl_ses=1:4%:6
        for freq=1:5
            tic
            dataCovs_op=squeeze(ave_hilcov_option3(:,subj,dl_ses,freq,:,:));
            [penalizationIn_op3(subj,dl_ses,freq),penalizationOut_op3(subj,dl_ses,freq),minDev_op3(subj,dl_ses,freq)]=...
            penaltyselection(SC,allLambdas,allLambdasOut,dataCovs_op);
            toc
        end
    end
end
% 28000/3600= 8 h x 12 = 72 hours
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
save('penaltyselection_op3_subj_1_dl_ses_1-4.mat','penalizationOut_op1','penalizationIn_op1','minDev_op1');

% fitprecision
% hilbert_dataCov_all=nan(12,2,12,5,894,894);
X_op3=nan(12,2,12,5,894,894);
for subj=2%:2
    for dl_ses=1%:6
        for freq=1:5
            penalizationIn=penalizationIn_op3(subj,dl_ses,freq);
            penalizationOut=penalizationOut_op3(subj,dl_ses,freq);
            for tr=1:12
                tic
                ses=1+2*(dl_ses-1);
                dataCov=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
                [X_op3(ses,subj,tr,freq,:,:)] = fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov);
                ses=2+2*(dl_ses-1);
                dataCov=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
                [X_op3(ses,subj,tr,freq,:,:)] = fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov);
                toc
            end
        end
    end
end
save('X_op3_subj_2_dl_ses_1.mat','X_op3','-v7.3');
% save('X_op3_subj_1_dl_ses_1.mat','X_op3');
% teest
% tic;xxx=fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov);toc; % 102 s
% estimated for all 100*2*5*6*2/3600 = 3 hours

%% option 4
% By comparison, if you just take 3 trials and 1 subject in 1 condition, 
% and do cross-validation what is the penalty inside and outside.
% penaltyselection for 1 subject

% In 1 subject first
% hilbert_dataCov_all=nan(12,2,12,5,894,894);
cd ../../Cleaned_data/hilbert_datacov
ave_hilcov_option4=nan(5,3,894,894); % 5 frequency x 3 ensamble x 894 x 894
ses=1; subj=1; tr=1:3;
for freq=1:5
    ave_hilcov_option4(freq,:,:,:)=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
end

penalizationIn_op4=nan(5,1);
penalizationOut_op4=nan(5,1);
minDev_op4=nan(5,1);
for freq=1:5
    tic
    dataCovs_op=squeeze(ave_hilcov_option4(freq,:,:,:));
    [penalizationIn_op4(freq),penalizationOut_op4(freq),minDev_op4(freq)]=...
    penaltyselection(SC,allLambdas,allLambdasOut,dataCovs_op);
    toc
end
save('penaltyselection_op4.mat','penalizationIn_op4','penalizationOut_op4','minDev_op4')

%% Compare subj=1, ses=1:2 of option 1 and 3
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov

condi=1;freq=5;
penalizationIn_op1(condi,freq)
penalizationOut_op1(condi,freq)
minDev_op1(condi,freq)

subj=1;dl_ses=1;freq=5;
penalizationIn_op3(subj,dl_ses,:)
penalizationOut_op3(subj,dl_ses,:)
minDev_op3(subj,dl_ses,:)


load('penaltyselection_op1.mat')
load('penaltyselection_op3.mat')

load('X_op1_ses_1-2.mat')
load('X_op3_subj_1_dl_ses_1.mat')

ses=1; subj=1;
freq=3;
tr=11;
% compare the precision output from op1 and op3
cmin=-0.001;cmax=0.001;
figure;
clf;
subplot(131);
imagesc(squeeze(X_op1(ses,subj,tr,freq,:,:)));colormap('jet');colorbar;
title(['option 1 ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr) ' freq ' num2str(freq)]);
clim([cmin cmax])

subplot(132);
imagesc(squeeze(X_op3(ses,subj,tr,freq,:,:)));colormap('jet');colorbar;
title(['option 3 ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr) ' freq ' num2str(freq)]);
clim([cmin cmax])

subplot(133);
imagesc(squeeze(X_op1(ses,subj,tr,freq,:,:))-squeeze(X_op3(ses,subj,tr,freq,:,:))); colormap('jet');colorbar;
title(['option 1-3 ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr) ' freq ' num2str(freq)]);
clim([cmin cmax])

% sum(diag(logical(squeeze(X_op3(ses,subj,tr,freq,:,:)))))
% (sum(logical(squeeze(X_op1(ses,subj,tr,freq,:,:))),'all')-894)/2
% sum(triu(logical(squeeze(X_op1(ses,subj,tr,freq,:,:)))),'all')-894

% compare 5 frequency in op1 and op3
figure;clf;
ses=1; subj=1;
tr=11;
cmin=-0.001;cmax=0.001;
for freq=1:5
    subplot(2,5,freq);
    imagesc(squeeze(X_op1(ses,subj,tr,freq,:,:)));colormap('jet');colorbar;
    title(['option 1 freq ' num2str(freq)]);
    clim([cmin cmax]);
    subplot(2,5,5+freq);
    imagesc(squeeze(X_op3(ses,subj,tr,freq,:,:)));colormap('jet');colorbar;
    title(['option 3 freq ' num2str(freq)]);
    clim([cmin cmax]);
end
sgtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr)]);

figure;clf;
ses=1; subj=2;
tr=11;
cmin=-0.001;cmax=0.001;
for freq=1:5
    subplot(2,5,freq);
    imagesc(squeeze(X_op1(ses,subj,tr,freq,:,:)));colormap('jet');colorbar;
    title(['option 1 freq ' num2str(freq)]);
    clim([cmin cmax]);
    subplot(2,5,5+freq);
    imagesc(squeeze(X_op3_ses_1_2(ses,subj,tr,freq,:,:)));colormap('jet');colorbar;
    title(['option 3 freq ' num2str(freq)]);
    clim([cmin cmax]);
end
sgtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr)]);

% compare the number of edges
addpath(genpath('/home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util'));
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('X_op1_ses_1-2.mat')
% load('X_op3_subj_1_dl_ses_1.mat')
% X_op3_ses_1_2=nan(12,2,12,5,894,894);
% X_op3_ses_1_2(1:2,1,:,:,:,:)=X_op3(1:2,1,:,:,:,:);
% load('X_op3_subj_2_dl_ses_1.mat');
% X_op3_ses_1_2(1:2,2,:,:,:,:)=X_op3(1:2,2,:,:,:,:);
% save('X_op3_ses_1_2.mat','X_op3_ses_1_2','-v7.3');
load('X_op3_ses_1_2.mat')

NedgeIn_op1=nan(12,2,12,5);
NedgeOut_op1=nan(12,2,12,5);
NedgeIn_op3=nan(12,2,12,5);
NedgeOut_op3=nan(12,2,12,5);
for ses =1:2
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                newG = reduce2nNetwork(logical(squeeze(X_op1(ses,subj,tr,freq,:,:))));
                X_reduced_op1(ses,subj,tr,freq,:,:)=newG;
                NedgeIn_op1(ses,subj,tr,freq) = sum(sum(newG.*triu(SC,1)));
                NedgeOut_op1(ses,subj,tr,freq) =  sum(sum(newG.*triu(~SC,1)));
                newG = reduce2nNetwork(logical(squeeze(X_op3_ses_1_2(ses,subj,tr,freq,:,:))));
                X_reduced_op1(ses,subj,tr,freq,:,:)=newG;
                NedgeIn_op3(ses,subj,tr,freq) = sum(sum(newG.*triu(SC,1)));
                NedgeOut_op3(ses,subj,tr,freq) =  sum(sum(newG.*triu(~SC,1)));
            end
        end
    end
end

% fraction
figure;
clf;
for ses=1:2
    for subj=1:2
        for tr=1:12
            if subj==1;
                subplot(4,12,(ses-1)*12+tr)
            elseif subj==2;
                subplot(4,12,(ses+1)*12+tr)
            end
            hold on;
            plot(1:5,squeeze(NedgeIn_op1(ses,subj,tr,:))/4985,'g');
            plot(1:5,squeeze(NedgeOut_op1(ses,subj,tr,:))/(99681-4985),'r');
            plot(1:5,squeeze(NedgeIn_op3(ses,subj,tr,:))/4985,'b');
            plot(1:5,squeeze(NedgeOut_op3(ses,subj,tr,:))/(99681-4985),'m');
            ylabel('n edges (prt)')
            xlabel('frequency band')
            title(['ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr)]);
            xticks([1:5])
            xticklabels(bandlabels)
            ylim([0 0.8])
            hold off;
            if ses==1 && subj==1 && tr==1;
                legend({'op1 in','op1 out','op3 in', 'op3 out'},'location','northwest');
            end
        end
    end
end



% compare op1 (across subj) and op3 (across condition)
std_subj_NedgeIn_op1=nan(2,5);
std_subj_NedgeOut_op1=nan(2,5);
mean_subj_NedgeIn_op1=nan(2,5);
mean_subj_NedgeOut_op1=nan(2,5);
for subj=1:2
    for freq=1:5
        % std(reshape(NedgeIn_op1(1:2,subj,:,freq),[],1))
        std_subj_NedgeIn_op1(subj,freq)=std(NedgeIn_op1(1:2,subj,:,freq),0,'all');
        std_subj_NedgeOut_op1(subj,freq)=std(NedgeOut_op1(1:2,subj,:,freq),0,'all');
        mean_subj_NedgeIn_op1(subj,freq)=mean(NedgeIn_op1(1:2,subj,:,freq),'all');
        mean_subj_NedgeOut_op1(subj,freq)=mean(NedgeOut_op1(1:2,subj,:,freq),'all');
    end
end


std_condi_NedgeIn_op3=nan(4,5);
std_condi_NedgeOut_op3=nan(4,5);
mean_condi_NedgeIn_op3=nan(4,5);
mean_condi_NedgeOut_op3=nan(4,5);
for condi=1:4
    for freq=1:5
        std_condi_NedgeIn_op3(condi,freq)=std(NedgeIn_op3(1:2,:,[1:3]+3*(condi-1),freq),0,'all');
        std_condi_NedgeOut_op3(condi,freq)=std(NedgeOut_op3(1:2,:,[1:3]+3*(condi-1),freq),0,'all');
        mean_condi_NedgeIn_op3(condi,freq)=mean(NedgeIn_op3(1:2,:,[1:3]+3*(condi-1),freq),'all');
        mean_condi_NedgeOut_op3(condi,freq)=mean(NedgeOut_op3(1:2,:,[1:3]+3*(condi-1),freq),'all');
    end
end

% try X*dataCov and see if = identity
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('hilbert_dataCov_all.mat')
load('X_op3_subj_1_dl_ses_1.mat')
subj=1;dl_ses=1;freq=1;
tr=1;
ses=1+2*(dl_ses-1);
dataCov=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
X=squeeze(X_op3(ses,subj,tr,freq,:,:));
% imagesc(X*dataCov);colorbar;
imagesc(X.*dataCov);colorbar;
% imagesc(X);colorbar;
% imagesc(dataCov);colorbar;
% imagesc(inv(dataCov));colorbar;
% imagesc(eye(894));colorbar;
% A=randi(10,10);
% imagesc(A);colorbar;
% imagesc(inv(A));colorbar;
% imagesc(A*inv(A));colorbar;
% imagesc(A*A');colorbar;
vlim=1;
clim([-1*vlim vlim]);
title('X.*dataCov')
subtitle(['Option 3: ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr)'])
% try real2Complex
covMat_complex = real2Complex(X, 0);

mean(diag(X.*dataCov))

network = ggmFitHtf(dataCov+ eye(length(dataCov)) *min(allLambdas)  * max(max(triu(abs(dataCov),1))),GforFit_new);

open inv


%% test X * cov and ggmFitHtf
cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/AGL
open estBestPenalizationQUIC.m

%%
% /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util/genCov_CMVN_SC.m
% /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util/normalizeCSD.m
% /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util/real2Complex.m


%% Other
% two papers
% drawmesh display the motor connections