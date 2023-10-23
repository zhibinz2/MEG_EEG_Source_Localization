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
% passbands = [1 3; 3.5 6.5; 7 10; 10.5 13.5; 14 20; 21 29; 30 49.5];
% bandlabels = {'Delta', 'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2', 'Gamma'};
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

for ses=1:12
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

                % ampdata = abs(hilbertdata); 
                % angledata = angle(hilbertdata);

                sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));
                sourceDataReal = sourceDataReal*(1/mean(abs(sourceDataReal(:)))); % normalize data

                hilbert_dataCov{freq} = cov(sourceDataReal');

                save(['./hilbert_datacov/' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat'],'hilbert_dataCov');
            end
        end
    end
    toc
end
 
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
figure;imagesc(SC);colorbar;colormap('jet')

% GforFit =[double(SC),double(SC) ; double(SC), double(SC)]; % boolean
% figure;imagesc(GforFit);colorbar;colormap('jet')
% 
% n_ensam=size(data,1); % number of ensambles


%% Organize hilber_dataCov_all
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

hilbert_dataCov_all=nan(12,2,12,5,894,894);
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
% save('hilbert_dataCov_all.mat','hilbert_dataCov_all');

%% option 1 By condition
% For any one frequency and condition, there are 72 covariance.  
% They come from  24 sessions (12 subjects x 2).   
% Lets break this up into 3 groups, each of which has 24 covariance matrix, drawn 1 from each subject and session. 
% You can think of it as grabbing the first trial for all sessions and subjects and making them 1 group, 
% the 2nd one is the 2nd group, 3rd one os the 3rd group.
% If we average those 24 together, then for each frequency band and condition we have 3 covariance matrix. 
% - put these 3 matrix into penalty selection and do 3 x cross validation to select penalty for that frequency band x condition.
% If we do it this way, we will have to depend on evaluating weights to compare conditions.

% Option 1 requires doing 4 x 5 cross-validation 
ave_hilcov_option1=nan(3,4,5,894,894); % 3 ensamble x 4 condition x 5 freq
for ensam=1:3
    for freq=1:5
        tic
        condi=1; % uncouple
        ave_hilcov_option1(ensam,condi,freq,:,:)=squeeze(mean(squeeze(mean( ...
            hilbert_dataCov_all(:,:,ensam+0*3,freq,:,:),1)),1));
        condi=2; % leading
        ave_hilcov_option1(ensam,condi,freq,:,:)=squeeze(mean(cat(1, ...
            squeeze(hilbert_dataCov_all(:,1,ensam+3,freq,:,:)),...
            squeeze(hilbert_dataCov_all(:,2,ensam+2*3,freq,:,:))),1));
        condi=3; % following
        ave_hilcov_option1(ensam,condi,freq,:,:)=squeeze(mean(cat(1, ...
            squeeze(hilbert_dataCov_all(:,2,ensam+3,freq,:,:)),...
            squeeze(hilbert_dataCov_all(:,1,ensam+2*3,freq,:,:))),1));
        condi=4; % mutual
        ave_hilcov_option1(ensam,condi,freq,:,:)=squeeze(mean(squeeze(mean( ...
            hilbert_dataCov_all(:,:,ensam+3*3,freq,:,:),1)),1));
        toc
    end
end
% 15 s
size(ave_hilcov_option1)

addpath(genpath('/../../AdaptiveGraphicalLassoforParCoh/Simulations/util'))
addpath(genpath('../../AdaptiveGraphicalLassoforParCoh/AGL'))

% penaltyselection
penalizationIn_op1=nan(4,5);
penalizationOut_op1=nan(4,5);
minDev_op1=nan(4,5);
for condi=1:4
    for freq=1:5
        tic
        dataCovs_op=squeeze(ave_hilcov_option1(:,condi,freq,:,:));
        [penalizationIn_op1(condi,freq),penalizationOut_op1(condi,freq),minDev_op1(condi,freq)]=...
            penaltyselection(SC,allLambdas,allLambdasOut,dataCovs_op);
        toc
    end
end
% 13 hours
save('penaltyselection_op1.mat','penalizationOut_op1','penalizationIn_op1','minDev_op1');

% fitprecision
% hilbert_dataCov_all=nan(12,2,12,5,894,894);
X_op1=nan(12,2,12,5,894,894);
parfor ses=1:12
    for subj=1:2
        for tr=1:12
            tic
            for freq=1:5
                if ismember(tr,[1 2 3]) % uncouple
                    penalizationIn=penalizationIn_op1(1,freq);
                    penalizationOut=penalizationOut_op1(1,freq);
                    dataCov=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
                    [X_op1(ses,subj,tr,freq,:,:)] = fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov);
                elseif (ismember(tr,[4:6]) & subj==1) | (ismember(tr,[7:9]) & subj==2) % leading
                    penalizationIn=penalizationIn_op1(2,freq);
                    penalizationOut=penalizationOut_op1(2,freq);
                    dataCov=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
                    [X_op1(ses,subj,tr,freq,:,:)] = fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov);
                elseif (ismember(tr,[4:6]) & subj==2) | (ismember(tr,[7:9]) & subj==1) % following
                    penalizationIn=penalizationIn_op1(3,freq);
                    penalizationOut=penalizationOut_op1(3,freq);
                    dataCov=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
                    [X_op1(ses,subj,tr,freq,:,:)] = fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov);
                else ismember(tr,[10:12])
                    penalizationIn=penalizationIn_op1(4,freq);
                    penalizationOut=penalizationOut_op1(4,freq);
                    dataCov=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
                    [X_op1(ses,subj,tr,freq,:,:)] = fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov);
                end
            end
            toc
        end
    end
end
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
save('X_op1.mat','X_op1')



%% Option 2:
% On each session there are 3 trials in 1 condition.  
% Get 3 covariance matrix, put them in cross-validation to select penalty.  
% Run all 3 trials at that penalty.  This would require doing 72 x 5 cross-validation models.  
% But, we would be less dependent on weights.
% Option 2 requires 2 x 12 x 4 x 5 cross-validation (not 72 x 5, estimated 33-38 days).
% I would vote for Option 1.

%% option 3 By subject
% We should also save the deviance values for the final fit, because they tell us how good a model it is.
% Option 3: Choose the penalties for each subject in each frequency band for all conditions.
% Each subject produces 24 trials, 3 in each condition.  So take 1 trial of each condition (8 trials 4 from each day) 
% and average covariance.  now you have 3 covariance matrices.  run the penalty selection function. .
% Go back and fit  each trial for that subject using that penalty.
% This way you would have to 12 x 5 penality optimizations.  But, I think this may be better, because the subject to subject differences may be greater than the conditon differences.
% Can you try this for the same frequency you are testing for option 1 in 1 subject.

% In 1 subject first
% hilbert_dataCov_all=nan(12,2,12,5,894,894);
cd ../../Cleaned_data/hilbert_datacov
ave_hilcov_option3=nan(3,2,6,5,894,894); % 3 ensamble x 2 subjects x 6 double-sessions x 5 frequency x 894 x 894
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
penalizationIn_op1=nan(2,6,5);
penalizationOut_op1=nan(2,6,5);
minDev_op1=nan(2,6,5);
parfor subj=1%:2
    for dl_ses=1%:6
        for freq=1:5
            tic
            dataCovs_op=squeeze(ave_hilcov_option3(:,subj,dl_ses,freq,:,:));
            [penalizationIn_op1(subj,dl_ses,freq),penalizationOut_op1(subj,dl_ses,freq),minDev_op1(subj,dl_ses,freq)]=...
            penaltyselection(SC,allLambdas,allLambdasOut,dataCovs_op);
            toc
        end
    end
end




%% test X * cov and ggmFitHtf
cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/AGL
open estBestPenalizationQUIC.m

%% Other
% two papers
% drawmesh display the motor connections