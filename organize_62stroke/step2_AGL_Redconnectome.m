%% downsample filter hilbert -> coh and partial coh
% refer to 
% open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/hilbert2cov.m
clear
Fs=1000;
srnew = 200;
downsample = Fs/srnew;

passbands = [1 3; 3.5 6.5; 7 10; 10.5 13.5; 14 20; 21 29];
bandlabels = {'Delta','Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};

attenuation=60;
filt_ds=cell(1,length(passbands));
for freqBand=1:6
    % Select frequency
    passFreq1 = passbands(freqBand,1);
    passFreq2 = passbands(freqBand,2);
    filt_d = designfilt('bandpassiir','FilterOrder',20, ...
    'PassbandFrequency1',passFreq1,'PassbandFrequency2',passFreq2, ...
    'StopbandAttenuation1',attenuation,'PassbandRipple',0.2, ...
    'StopbandAttenuation2',attenuation,'SampleRate',srnew);
    filt_ds{freqBand}=filt_d;
end

% AGL variables
allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
n_Lambdas=length(allLambdas); % number of lambda values
min_LamdaIn = min(allLambdas);

cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source
load('corti_ave_source_labl.mat'); % every trial is the same, load one trial

% penaltyselection for 1 subject
addpath(genpath('/../../AdaptiveGraphicalLassoforParCoh/Simulations/util'))
addpath(genpath('../../AdaptiveGraphicalLassoforParCoh/AGL'))
cd ../../Virtual-Tractography/ForZhibin/processed_data/
load('scale250_Connectome.mat','fc');
sum(triu(logical(fc),1),"all") % 6386 edges
SC=logical(fc(corti_ave_source_labl,corti_ave_source_labl));
sum(triu(SC,1),"all") % 4990 edges

% compare full connectome to fc
figure;imagesc(logical(fc));
cd /home/zhibinz2/Documents/GitHub/STROKE_P61/getLesionMaskConnectome_nonflip_p61
load('Connectome_p61.mat', 'fullConnectome_p61')
figure;imagesc(logical(squeeze(fullConnectome_p61(61,:,:)))); % look similar to fc though
sum(triu(logical(squeeze(fullConnectome_p61(61,:,:))),1),"all") % 5814 edges

% extract indivisual reduced sc
load('Connectome_p61.mat', 'redConnectome_p61')