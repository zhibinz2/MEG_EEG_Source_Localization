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
% cd ../../Virtual-Tractography/ForZhibin/processed_data/
% load('scale250_Connectome.mat','fc');
% sum(triu(logical(fc),1),"all") % 6386 edges
% SC=logical(fc(corti_ave_source_labl,corti_ave_source_labl));
% sum(triu(SC,1),"all") % 4990 edges

% % compare full connectome to fc
% figure;imagesc(logical(fc));
% cd /home/zhibinz2/Documents/GitHub/STROKE_P61/getLesionMaskConnectome_nonflip_p61
% load('Connectome_p61.mat', 'fullConnectome_p61')
% figure;imagesc(logical(squeeze(fullConnectome_p61(61,:,:)))); % look similar to fc though
% sum(triu(logical(squeeze(fullConnectome_p61(61,:,:))),1),"all") % 5814 edges

% extract indivisual reduced sc
cd /home/zhibinz2/Documents/GitHub/STROKE_P61/getLesionMaskConnectome_nonflip_p61
load('Connectome_p61.mat', 'redConnectome_p61')
SC_p61=nan(61,448,448);
for p=1:61
    redC_mat=squeeze(redConnectome_p61(p,:,:));
    SC_p61(p,:,:)=logical(redC_mat(corti_ave_source_labl,corti_ave_source_labl));
end
% imagesc(logical(redC_mat(corti_ave_source_labl,corti_ave_source_labl)));colorbar();
% sum(triu(logical(redC_mat(corti_ave_source_labl,corti_ave_source_labl)),1),"all") % 5790 edges

%% penalty selection and fit precision
subj_files=[0:1:60];
for f=29:35
    tic
    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source
    load([num2str(subj_files(f)) '.mat'],'corti_source_data','subject_ID', ...
        'chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');
    display(['start processing subject file: ' num2str(subj_files(f)) '.mat']);

    downsample_data=resample(double(corti_source_data),1,downsample,'Dimension',1);
    
    datapermuted=cell(1,6);
    stroke_Cov=cell(1,6); % coh
    stroke_Pcov=cell(1,6); % partial coh
    penalizationIn_op=nan(1,6);
    penalizationOut_op=nan(1,6);
    minDev_op=nan(1,6);
    for freq=1:6
        filtered_data = filter(filt_ds{freq},downsample_data);
        hilbertdata = hilbert(filtered_data');
        sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));
        sourceDataReal = [sourceDataReal*(1/mean(abs(sourceDataReal(:))))]'; % normalize data
        % split into 4 ensambles: 4 x #source x #(samples/4)
        n_split=2; n_sr=size(sourceDataReal,2);
        sam_len=size(sourceDataReal,1);
        sam_size=floor(sam_len/n_split); sam_range=1:sam_size;
        sourceDataReal = sourceDataReal(1:n_split*sam_size,:);
        datareshaped = reshape(sourceDataReal, sam_size, n_split, n_sr);
        datapermuted{freq} = permute(datareshaped,[2,3,1]); % split into 2 ensambles: 2x896x9300
        % compute covariance
        stroke_Cov{freq} = cov(sourceDataReal);
        for n=1:n_split
            datapermuted_cov{freq}(n,:,:) = cov([squeeze(datapermuted{freq}(n,:,:))]');
        end

        clear corti_source_data filtered_data  hilbertdata
        clear sourceDataReal datareshaped 
        clear SC

        % Individual SC from redConnectome
        SC=squeeze(SC_p61(f,:,:));

        % AGL
        cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
        dataCovs_op=squeeze(datapermuted_cov{freq});
        [penalizationIn_op(freq),penalizationOut_op(freq),minDev_op(freq)]=penaltyselection( ...
        SC,allLambdas,allLambdasOut,dataCovs_op);
        [stroke_Pcov{freq}] = fitprecision( ...
        SC,penalizationIn_op(freq),penalizationOut_op(freq),min_LamdaIn,stroke_Cov{freq});
        
        clear datapermuted_cov dataCovs_op 

    end
    clear downsample_data datapermuted

    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_61_source_hilbert_redC
    save([num2str(subject_ID) '.mat'],'stroke_Cov', 'stroke_Pcov', ...
        'penalizationIn_op','penalizationOut_op','minDev_op', ...
        'subject_ID','chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');

    display(['Complete one file: ' num2str(subj_files(f)) '.mat ********************'])
    toc
end