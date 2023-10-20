%% get the common variables 
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

load('source_rr.mat');
load('Lausanne2008_fsaverageDSsurf_60_125_250.mat')
% label the sources
Vertex=Brain.Vertex;
load('parcels.mat') % This is the labels
% Anni's labeling method
x_shift=(max(Vertex(:,1))-max(source_rr(:,1))*1e3)/2+(min(Vertex(:,1))-min(source_rr(:,1))*1e3)/2;
y_shift=(max(Vertex(:,2))-max(source_rr(:,2))*1e3)/2+(min(Vertex(:,2))-min(source_rr(:,2))*1e3)/2;
z_shift=(max(Vertex(:,3))-max(source_rr(:,3))*1e3)/2+(min(Vertex(:,3))-min(source_rr(:,3))*1e3)/2;
source_x=source_rr(:,1) * 1e3 + x_shift;
source_y=source_rr(:,2) * 1e3 + y_shift;
source_z=source_rr(:,3) * 1e3 + z_shift;
source_xyz=[source_x source_y source_z];
num_source=size(source_xyz,1);
source_fsaverage = source_xyz+127.5; % 127.5 is based on the fsaverage volume being 256 x 256 x 256
source_labels=zeros(num_source,1);
for i = 1:length(source_fsaverage)
    vox = ceil(source_fsaverage(i,:));
    inds              = sub2ind([size(parcels)], vox(1), vox(2), vox(3));
    label             = parcels(inds); 
    source_labels(i) = label;
end

% load forward matrix
load('leadfield.mat');

% inverse model
addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat, stat, reconstructed] = inversemodel(leadfield,'prctile',1);
% leakage=inversemat*leadfield;

%% loop through all data to compute correlation of the recon EEG with original EEG 
data_path= '../../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 
tic
corr_32_all=zeros(numSes,2,12,32);
for r=1:numSes
    clear runid
    runid = num2str(seeds(r,:));
    % data=load([data_path  'clean_' runid '.mat']);
    data=load([data_path  'clean_' runid '.mat'],'dataL','dataR');
    for subj=1:2
        if subj==1
            preprocessed_eeg=data.dataL;
        else
            preprocessed_eeg=data.dataR;
        end
        for tr=1:12
            EEG_ori=preprocessed_eeg{tr}(:,1:32)';
            source_data=inversemat*EEG_ori;
            EEG_recon=leadfield*source_data;

            for ch=1:32
                R=corrcoef(EEG_recon(ch,:)',EEG_ori(ch,:)');
                corr_32_all(r,subj,tr,ch)=R(1,2);
            end
        end
    end
end
toc
% 440 s 

mean(corr_32_all,4)


%% loop throught all source data to compute pca and aggregate
for r=1:numSes
    clear runid data preprocessed_eeg EEG_ori source_data EEG_recon leakage
    runid = num2str(seeds(r,:));
    % data=load([data_path  'clean_' runid '.mat']);
    data=load([data_path  'clean_' runid '.mat'],'dataL','dataR');
    for subj=1:2
        if subj==1
            preprocessed_eeg=data.dataL;
        else
            preprocessed_eeg=data.dataR;
        end
        for tr=1:12
            clear EEG_ori source_data EEG_recon leakage
            tic
            EEG_ori=preprocessed_eeg{tr}(:,1:32)';
            source_data=inversemat*EEG_ori;
            EEG_recon=leadfield*source_data;

            fra_eigenvalues=zeros(1,max(unique(source_labels)));
            agr_source_data=[];
            ave_source_coor=[];
            ave_source_label=[];
            for sr=1:max(unique(source_labels))
                I=find(source_labels==sr);
                if ~isempty(I)
                    [COEFF, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
                    fra_eigenvalues(sr)=LATENT(1)/sum(LATENT);
                    if fra_eigenvalues(sr) > 0.5
                        agr_source_data=[agr_source_data SCORE(:,1)];
                        ave_source_coor=[ave_source_coor; mean(source_fsaverage(I,:),1)];
                        ave_source_label=[ave_source_label; sr];
                    end
                end
            end
            save([data_path 'source_data/' num2str(runid) '/subj' num2str(subj) '_tr_' num2str(tr)  '.mat'], ...
                            'fra_eigenvalues','agr_source_data','ave_source_coor','ave_source_label');
            toc
        end
    end
end
% 9 h

%% try remove subcortical and "zeros marked" sources not mapped and saved the aggreated source data
sum(fra_eigenvalues==0)
sum(fra_eigenvalues>0.50)
roiNames_250(scale250_subcortROIs(2:end))
fra_eigenvalues(scale250_subcortROIs(2:end))=0
% 463-15=448 cortical sources

% for ses=1:12
%     mkdir(num2str(seeds(ses,:)));
% end


cd /home/zhibinz2/Documents/GitHub/Cleaned_data/cortical_source_data
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
load('/home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_32chan_pca/Lausanne2008_fsaverageDSsurf_60_125_250.mat')
corti_fra_eigenvalues=cell(12,2,12);
corti_ave_source_coor=cell(12,2,12);
corti_ave_source_labl=cell(12,2,12);
for ses=1:12
    for subj=1:2
        for tr=1:12
            tic
            clear i agr_source_data ave_source_coor ave_source_label
            load(['../source_data/' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat'])
            marked_fra_eigenvalues=fra_eigenvalues(ave_source_label);
            bool_temp=zeros(1,length(ave_source_label));
            for i=1:length(ave_source_label)
                clear tmp
                tmp=ave_source_label(i);
                if ismember(tmp, scale250_subcortROIs)
                    bool_temp(i)=1;
                end
            end
            clear ind
            ind=find(bool_temp);
            agr_source_data(:,ind)=[];
            marked_fra_eigenvalues(ind)=[];
            ave_source_coor(ind,:)=[];
            ave_source_label(ind)=[];
            corti_fra_eigenvalues{ses,subj,tr}=marked_fra_eigenvalues;
            corti_ave_source_coor{ses,subj,tr}=ave_source_coor;
            corti_ave_source_labl{ses,subj,tr}=ave_source_label;
            save(['./' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat'],'agr_source_data')
            toc
        end
    end
end
save('corti_fra_eigenvalues.mat','corti_fra_eigenvalues');
save('corti_ave_source_coor.mat','corti_ave_source_coor');
save('corti_ave_source_labl.mat','corti_ave_source_labl');


%% count number of souces > 0.5 remain
num_source_all=nan(12,2,12);
for ses=1:12
    for subj=1:2
        for tr=1:12
            num_source_all(ses,subj,tr)=length(corti_ave_source_labl{ses,subj,tr});
        end
    end
end
num_source_all

min(num_source_all,[],'all') % 448-430 less than 18 ROIs has fraction of 1st eigenvalue < 0.5

max(num_source_all,[],'all')



%% loop through all data to compute fraction of the first eigenvalue
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/source_data
% ses=1; subj=1; tr=5;
fra_eigenvalues_all=zeros(12,2,12,463);
for ses =1:12
    for subj =1:2
        tic
        for tr=1:12
            clear fra_eigenvalues
            load(['./' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat'])
            fra_eigenvalues_all(ses,subj,tr,:)=fra_eigenvalues;
%             prctile(fra_eigenvalues, 25)
%             if median(fra_eigenvalues) < 0.7
%                 display(['ses ' num2str(ses) ' subj ' num2str(subj) ' tr ' num2str(tr) '<0.7'])
%             end
%             plot(1:max(unique(source_labels)), fra_eigenvalues,'.','MarkerSize',12);
%             yline(0.5);yline(0.7);
%             xlabel('ROI labels');ylabel('fraction of 1st eigenvalue')
        end
        toc
    end
end
% 0.5 h?

fra_median=zeros(12,2,12); % mean 0.798
fra_25prct=zeros(12,2,12); % mean 0.6859
for ses =1:12
    for subj =1:2
        for tr=1:12
            fra_median(ses,subj,tr)=median(fra_eigenvalues_all(ses,subj,tr,:));
            fra_25prct(ses,subj,tr)=prctile(fra_eigenvalues_all(ses,subj,tr,:), 25);
        end
    end
end

figure;
clf;
legends=cell(1,2);
plot(1:288, reshape(fra_median,[],1),'r.','MarkerSize',10);
legends{1}='median';
hold on;
plot(1:288, reshape(fra_25prct,[],1),'b.','MarkerSize',10);
legends{2}='25 percentile';
ylabel('fraction of 1st eigenvalue');xlabel('all 288 trials');
legend(legends);

%% agr_source_data - hilbert - coh
% https://www.mathworks.com/help/signal/ref/ellip.html

% example of bandpass filter
sampl_rate=2000;
[A,B,C,D] = ellip(20,0.25,60,[14 20]/(sampl_rate/2));
d = designfilt('bandpassiir','FilterOrder',20, ...
    'PassbandFrequency1',14,'PassbandFrequency2',20, ...
    'StopbandAttenuation1',60,'PassbandRipple',0.25, ...
    'StopbandAttenuation2',60,'SampleRate',2000);

sos = ss2sos(A,B,C,D);
fvt = fvtool(sos,d,'Fs',sampl_rate);% visualize 
legend(fvt,'ellip','designfilt') 

dataIn = randn(1000,1);plot(1:1000,dataIn,'g');hold on;
dataOut = filter(d,dataIn); plot(1:1000,dataOut,'r');hold off;


% example of a bandstop filter
[b,a] = ellip(3,5,50,[0.01 0.5],'stop');
freqz(b,a)
clf;
dataIn = randn(1000,1);plot(1:1000,dataIn,'g');hold on;
dataOut = filter(b,a,dataIn); plot(1:1000,dataOut,'r');hold off;


% ChatGPT convert Python2matlab
sampl_rate=2000;
srnew = 200;
downsample = 10;
passbands = [1 3; 3.5 6.5; 7 10; 10.5 13.5; 14 20; 21 29; 30 49.5];
% stopbands = [0.75 3.25; 3.25 6.75; 6.75 10.25; 10.25 13.75; 13.75 20.5; 20.5 29.5; 29.5 50];
bandlabels = {'Delta', 'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2', 'Gamma'};
% My adapted code  
freqBand=2;
attenuation=60;
passFreq1 = passbands(freqBand,1);
passFreq2 = passbands(freqBand,2);
% [A,B,C,D] = ellip(20,0.25,attenuation,passbands(freqBand,:)/(srnew/2));
d = designfilt('bandpassiir','FilterOrder',20, ...
    'PassbandFrequency1',passFreq1,'PassbandFrequency2',passFreq2, ...
    'StopbandAttenuation1',attenuation,'PassbandRipple',0.2, ...
    'StopbandAttenuation2',attenuation,'SampleRate',srnew);
% sos = ss2sos(A,B,C,D);
% fvt = fvtool(sos,d,'Fs',srnew);
% legend(fvt,'ellip','designfilt')

% dataIn = random_sr_data(151000:153000);clf;plot(1:2001,dataIn,'g');hold on;
% dataOut = filter(d,dataIn); plot(1:2001,dataOut,'r');hold off;

% % test on single channel
% random_sr_data=agr_source_data(:,250);
% % random_sr_data=EEG_ori(20,:)';
% dataIn = random_sr_data(151001:153000); %clf;plot(1:2000,dataIn,'g');hold on;
% downsample_data=resample(double(dataIn),1,downsample,'Dimension',1); % plot(1:downsample:2000,downsample_data,'b');
% filterd_data = filter(d,downsample_data); clf; plot(1:downsample:2000,filterd_data,'r'); hold on;
% hilbertdata = hilbert(filterd_data'); 
% ampdata = abs(hilbertdata); plot(1:downsample:2000,ampdata,'m');
% legend(['original source data'],['downsampled'],['filtered'],['hilbert'])
% angledata = angle(hilbertdata);
% 
% % coherence / cross spectra
% coh=cov(hilbertdata');
% % correlation
% cor=corrcoef(ampdata')


% test on all channel
random_sr_data=agr_source_data;
dataIn = double(random_sr_data(151001:153000,:));
downsample_data=resample(double(dataIn),1,downsample,'Dimension',1);
filterd_data = filter(d,downsample_data);
hilbertdata = hilbert(filterd_data'); 
ampdata = abs(hilbertdata); 
% angledata = angle(hilbertdata);
% coherence / cross spectra
coh=cov(hilbertdata');
covMat = real2Complex(coh,1);% it doesn't convert to real 
% correlation
cor=corrcoef(ampdata');
figure;
clf;
imagesc(cor);colorbar;colormap("jet");
title(['Amplitude correlation of the Hilbert data from PCA source time series']);
subtitle(['Freq band - ' bandlabels{freqBand}]);

%% AGL pcoh
% load source data
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/cortical_source_data/20220713
load('subj1_tr_1.mat')
cd ../
load('corti_ave_source_labl.mat')
load('corti_ave_source_coor.mat')
source_labels=corti_ave_source_labl{1,1,1};
source_coor=corti_ave_source_coor{1,1,1};
samps=size(hilbertdata,2);
sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));    
sourceDataReal = sourceDataReal*(1/mean(abs(sourceDataReal(:)))); % normalize data
datareshaped = permute(reshape(sourceDataReal', 4,samps/4, size(sourceDataReal,1)),[1,3,2]);
% load SC
load('../../Virtual-Tractography/ForZhibin/processed_data/scale250_Connectome.mat');
SC_tr=fc(source_labels,source_labels);
figure;imagesc(fc);colorbar;colormap("jet");clim([0 0.2]);
figure;imagesc(logical(SC_tr));colorbar;colormap("jet");
% AGL
allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]); 
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
addpath(genpath('/../../AdaptiveGraphicalLassoforParCoh/Simulations/util'))
addpath(genpath('../../AdaptiveGraphicalLassoforParCoh/AGL'))
tic
clear networkPrecCompTrue penInCompTrue penOutCompTrue allDevsReturnTrue
[networkPrecCompTrue, penInCompTrue, penOutCompTrue,~,allDevsReturnTrue] ...
        = estBestPenalizationQUIC(datareshaped, logical(SC_tr), allLambdas,allLambdasOut, 0);
toc % 