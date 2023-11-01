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
addpath(genpath('/../../AdaptiveGraphicalLassoforParCoh/Simulations/util'))
addpath(genpath('../../AdaptiveGraphicalLassoforParCoh/AGL'))
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
save('penaltyselection_op3_subj_1_dl_ses_1-4.mat','penalizationOut_op3','penalizationIn_op3','minDev_op3');

% combine hnlb, hnba, hnl1
penalizationIn_op3_all=nan(2,6,5);
penalizationOut_op3_all=nan(2,6,5);
minDev_op3_all=nan(2,6,5);
load('penaltyselection_op3_subj_1_dl_ses_1-4.mat')
penalizationIn_op3_all(1,1:4,:)=penalizationIn_op3(1,1:4,:);
penalizationOut_op3_all(1,1:4,:)=penalizationOut_op3(1,1:4,:);
minDev_op3_all(1,1:4,:)=minDev_op3(1,1:4,:);
load('penaltyselection_op3_subj_1_dl_ses_5-6.mat')
penalizationIn_op3_all(1,5:6,:)=penalizationIn_op3(1,5:6,:);
penalizationOut_op3_all(1,5:6,:)=penalizationOut_op3(1,5:6,:);
minDev_op3_all(1,5:6,:)=minDev_op3(1,5:6,:);
load('penaltyselection_op3_subj_2_dl_ses_1-6.mat')
penalizationIn_op3_all(2,1:6,:)=penalizationIn_op3(2,1:6,:);
penalizationOut_op3_all(2,1:6,:)=penalizationOut_op3(2,1:6,:);
minDev_op3_all(2,1:6,:)=minDev_op3(2,1:6,:);
save('penaltyselection_op3_all.mat','penalizationOut_op3_all','penalizationIn_op3_all','minDev_op3_all');
% subj=2;ses_dl=1;freq=1; 

% fitprecision
% hilbert_dataCov_all=nan(12,2,12,5,894,894);
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
X_op3=nan(12,2,12,5,896,896);
for subj=1%:2
    for dl_ses=1:4%:6
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
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
save('X_op3_subj_1_dl_ses_1-4.mat','X_op3','-v7.3');
% save('X_op3_subj_1_dl_ses_1.mat','X_op3');
% teest
% tic;xxx=fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov);toc; % 102 s
% estimated for all 100*2*5*6*2/3600 = 3 hours

% % combine hnlb, hnba, hnl1
X_op3_all=nan(12,2,12,5,896,896);
load('X_op3_subj_1_dl_ses_1-4.mat')
for subj=1%:2
    for dl_ses=1:4%:6
        for freq=1:5
            for tr=1:12
                ses=1+2*(dl_ses-1);
                X_op3_all(ses,subj,tr,freq,:,:)=X_op3(ses,subj,tr,freq,:,:);
                ses=2+2*(dl_ses-1);
                X_op3_all(ses,subj,tr,freq,:,:)=X_op3(ses,subj,tr,freq,:,:);
            end
        end
    end
end
load('X_op3_subj_1_dl_ses_5-6.mat')
for subj=1%:2
    for dl_ses=5:6%:6
        for freq=1:5
            for tr=1:12
                ses=1+2*(dl_ses-1);
                X_op3_all(ses,subj,tr,freq,:,:)=X_op3(ses,subj,tr,freq,:,:);
                ses=2+2*(dl_ses-1);
                X_op3_all(ses,subj,tr,freq,:,:)=X_op3(ses,subj,tr,freq,:,:);
            end
        end
    end
end
load('X_op3_subj_2_dl_ses_1-6.mat')
for subj=2%:2
    for dl_ses=1:6%:6
        for freq=1:5
            for tr=1:12
                ses=1+2*(dl_ses-1);
                X_op3_all(ses,subj,tr,freq,:,:)=X_op3(ses,subj,tr,freq,:,:);
                ses=2+2*(dl_ses-1);
                X_op3_all(ses,subj,tr,freq,:,:)=X_op3(ses,subj,tr,freq,:,:);
            end
        end
    end
end
save('X_op3_all.mat','X_op3_all','-v7.3');

%% Examing number of edges
addpath('/home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util')
NedgeIn_op3=nan(12,2,12,5);
NedgeOut_op3=nan(12,2,12,5);
Reduced_X_op3_all=nan(12,2,12,5,448,448);
tic
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                X=squeeze(X_op3_all(ses,subj,tr,freq,:,:));
                newG = reduce2nNetwork(logical(X));
                NedgeIn_op3(ses,subj,tr,freq) = sum(sum(newG.*triu(SC,1)));
                NedgeOut_op3(ses,subj,tr,freq) =  sum(sum(newG.*triu(~SC,1)));
                Reduced_X_op3_all(ses,subj,tr,freq,:,:)=newG;
            end
        end
    end
end
toc % 29s
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
save('Reduced_X_op3_all.mat','Reduced_X_op3_all','-v7.3');


Complex_X_op3_all=nan(12,2,12,5,448,448);
tic
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                X=squeeze(X_op3_all(ses,subj,tr,freq,:,:));
                Complex_X_op3_all(ses,subj,tr,freq,:,:)=r2c(X);
            end
        end
    end
end
toc % 33s 
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
save('Complex_X_op3_all.mat','Complex_X_op3_all','-v7.3');


% examine fraction
for ses=1:12
    figure
    for subj=1:2
        for tr=1:12
            if subj==1;
                subplot(2,12,tr)
            elseif subj==2;
                subplot(2,12,12+tr)
            end
            hold on;
            plot(1:5,squeeze(NedgeIn_op3(ses,subj,tr,:))/4985,'b');
            plot(1:5,squeeze(NedgeOut_op3(ses,subj,tr,:))/(99681-4985),'m');
            ylabel('n edges (prt)')
            xlabel('frequency band')
            title(['ses ' num2str(ses) ' subj ' num2str(subj) ' trial ' num2str(tr)]);
            xticks([1:5])
            xticklabels(bandlabels)
            ylim([0 0.8])
            hold off;
            if subj==1 && tr==1 && subj==1;
                legend({'op3 in', 'op3 out'},'location','northwest');
            end
        end
    end
end


%% Examing subj=2;ses_dl=1;freq=1; (test)
% try X*dataCov and see if = identity
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('hilbert_dataCov_all.mat')
load('X_op3_all.mat')
subj=2;ses_dl=1;freq=1; 
tr=1;
ses=1+2*(ses_dl-1);
dataCov=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
X=squeeze(X_op3_all(ses,subj,tr,freq,:,:));
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
covMat_complex = real2Complex(X, 0); % no dimension reduction


mean(diag(X.*dataCov))

newG = reduce2nNetwork(logical(X));
GforFit_new = [newG, newG.*double(~eye(length(newG))); newG.*double(~eye(length(newG))), newG];
network = ggmFitHtf(dataCov+ eye(length(dataCov)) *min(allLambdas)  * max(max(triu(abs(dataCov),1))),GforFit_new);
% imagesc(network*dataCov);colorbar;
open inv

