% original data
cd ../../archive/EEG_StrokePatients_n61/

% combine files
cd /home/zhibinz2/Documents/GitHub/archive/EEG_StrokePatients_n61/HORR
load('HORR_DNDB78_20111214_14421.mat')
data1=Data;
load('HORR_DNDB78_20111214_14422.mat')
data2=Data;
load('HORR_DNDB78_20111214_14423.mat')
data3=Data;
Data=cat(2,data1,data2,data3);
save("HORR_DNDB78_20111214_all.mat","Data")

cd /home/zhibinz2/Documents/GitHub/archive/EEG_StrokePatients_n61/MAVP
load('MAVP_PWKZ25_20120118_14351.mat')
data1=Data;
load('MAVP_PWKZ25_20120118_14352.mat')
data2=Data;
load('MAVP_PWKZ25_20120118_14353.mat')
data3=Data;
Data=cat(2,data1,data2,data3);
save("MAVP_PWKZ25_20120118_all.mat",'Data')

cd /home/zhibinz2/Documents/GitHub/archive/EEG_StrokePatients_n61/MUNE
load('MUNE_XIUV73_20120215_12491.mat')
data1=Data;
load('MUNE_XIUV73_20120215_12492.mat')
data2=Data;
load('MUNE_XIUV73_20120215_12493.mat')
data3=Data;
Data=cat(2,data1,data2,data3);
save('MUNE_XIUV73_20120215_all.mat','Data')

cd /home/zhibinz2/Documents/GitHub/archive/EEG_StrokePatients_n61/SARB
load('SARB_QJQJ79_20120208_14231.mat')
data1=Data;
load('SARB_QJQJ79_20120208_14232.mat')
data2=Data;
load('SARB_QJQJ79_20120208_14233.mat')
data3=Data;
Data=cat(2,data1,data2,data3);
save('SARB_QJQJ79_20120208_all.mat',"Data")

cd /home/zhibinz2/Documents/GitHub/archive/EEG_StrokePatients_n61/STAP
load('STAP_EXNZ50_20111221_12381.mat')
data1=Data;
load('STAP_EXNZ50_20111221_12382.mat')
data2=Data;
load('STAP_EXNZ50_20111221_12383.mat')
data3=Data;
Data=cat(2,data1,data2,data3);
save('STAP_EXNZ50_20111221_all.mat','Data');



%%
% renamed accoording to post stroke days
open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/organize_62stroke/lesion_site_62subj_sort.xls
cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_reorganized
subject_ID=0;
subj_files=[0:1:60];
load([num2str(subject_ID) '.mat']);

% below refer to /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_spacing_ESCH/step1_EEG_prepocess.m
%% average re-reference
rData=Data-ones(size(Data,1),1)*mean(Data,1);

%% remove bad channels
% remove 49 peripheral channels
ch_peripheral=[241 242 243 238 239 240 ...
    244 245 246 247 251 256 91 102 111 120 133 145 165 174 187 199 208 216 229 233 237 236 235 234 ...
    232 228 217 209 200 188 175 166 156 146 134 121 112 103 92 82 255 250 ...
    257];
ch_peripheral=sort(ch_peripheral);
rData(ch_peripheral,:)=[];
rData=rData';
% detrend the EEG data (no padding needed)
detrend_data=detrend(rData,1);
% add paddings
padding=zeros(round(size(detrend_data,1)/10), size(detrend_data,2));
detrend_pad=cat(1,padding,detrend_data,padding);
% High pass
% addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab_zhibin/EEG/hnl'))
Fs=1000;
% high pass (no paddings needed)
Hd = makefilter(Fs,0.25,0.01,6,20,0); % for keeping readiness potential
pad_hp=filtfilthd(Hd,detrend_pad);
% Low pass (take a few seconds)
Hd = makefilter(Fs,50,51,6,20,0);  
pad_lp=filtfilthd(Hd,pad_hp);
% remove paddings
filtered_data=pad_lp((size(padding,1)+1):(size(padding,1)+size(detrend_data,1)),:);
% clear out memory before proceeding to next section
clearvars Data rData detrend_data detrend_pad pad_hp pad_lp
clearvars padding Hd
% organized into 1s epochs
nepochs=floor(length(filtered_data)/Fs);
nchans=size(filtered_data,2);
epochdata=zeros(nepochs,Fs,nchans);
for e = 1: nepochs
    epochdata(e,:,:)=filtered_data((1+(e-1)*Fs):(Fs+(e-1)*Fs),:);
end
clearvars e

%% identify bad chan and epochs by standard deviation
eegstd=squeeze(std(epochdata,[],1));
chanstd=sum(eegstd,2);
epochstd=sum(eegstd,1);

var_threshold = 2.5;  % normalized variance threshold to reject trials.
chan_threshold = 2.5;  % to identify % (smaller values are stricter)
corr_threshold = 0.4;  % threshold for identifying whether an ICA component contains eye movement. (smaller values are stricter)

% threshhold epochs by standard deviation criteria and remove them.
% in principle you could threshold channels this way too.  But, I
% think with 32 channels you need to avoid that.  With 128 you could.
badchan = find(epochstd / median(epochstd) > chan_threshold);
goodchan = setdiff(1:nchans,badchan);

badepoch = find(chanstd / median(chanstd) > var_threshold);
goodepoch = setdiff(1:nepochs, badepoch);

% remove bad epochs and connect back to timexchan 
good_epochdata=squeeze(epochdata(goodepoch,:,:));
good_filtered_data=zeros(length(goodepoch)*Fs,nchans);
for n = 1:length(goodepoch)
    good_filtered_data(((1+(n-1)*Fs):(Fs+(n-1)*Fs)),:)=squeeze(good_epochdata(n,:,:));
end

clearvars good_epochdata epochdata epochstd eegstd cahnstd
clearvars chan_threshold var_threshold filtered_data

%% extract good chanlocs
chanlocs = load('/home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_streamline_DAON/chanlocs.mat');
chanlocs = chanlocs.chanlocs;

% Add the original sequence index
for ch=1:257 % Channel 257 is Cz
    chanlocs(ch).urchan=ch;
end

% Create channel labels
ch_labels = cell(256,1);
for c = 1:256
    ch_labels{c}=num2str(c);
end
ch_labels{18}='Fp2'; ch_labels{37}='Fp1';
ch_labels{36}='F3'; ch_labels{224}='F4';
ch_labels{59}='C3'; ch_labels{183}='C4';
ch_labels{69}='T7'; ch_labels{202}='T8';
ch_labels{87}='P3'; ch_labels{153}='P4';
ch_labels{96}='P7'; ch_labels{170}='P8';
ch_labels{94}='LM'; ch_labels{190}='RM';
ch_labels{116}='O1'; ch_labels{150}='O2';
ch_labels{31}='NAS'; 
ch_labels{21}='Fz'; ch_labels{101}='Pz'; ch_labels{126}='Oz';

% remove 49 perpheral channels (then to check channel index for eyechans later)
central_chanlocs=chanlocs;
central_chanlocs(ch_peripheral)=[];

% remove Cz channel from ch_bad
ch_peripheral(end)=[];
ch_peripheral

% use "urchan" to get an array of the orginal index
urchan=extractfield(central_chanlocs,'urchan');

clearvars c ch chanloc

%% run ICA
% addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab_zhibin/EEG/hnl'))
tic
[icasig, A, W] = fastica(good_filtered_data'); 
toc % take a few minutes

% Assign dubious chans for correlation
eyechans = [248 252 253 67 249 254 73 230 226 225 219 231 227 218]; eyechans=sort(eyechans);
% the above eyechans  67    73   218   219   225   226   227   230   231   248   249   252   253   254
% they corresponds to 67    73   192   193   199   200   201   202   203   204   205   206   207   208 
% in the updated chanel index of goodchanlocs (icachans below)

% identify the updated indices for eyechan index 
icachans=nan(1,length(eyechans));
for c=1:length(eyechans)
    icachans(c)=find(urchan==eyechans(c));
end
% combine into dubious chans with badchan for ica
dubious_chans =unique([badchan icachans]);

% compute correlation with dubious chans (a few seconds)
corrs=zeros(length(dubious_chans),size(icasig,1));
for du = 1:length(dubious_chans)
    corrs(du, :)=corr(icasig',good_filtered_data(:,dubious_chans(du)));
end
corrsmax = max(corrs);

% figure;
% subplot(211);imagesc(corrs);colormap('jet');xlabel('components');ylabel('dubious chans'); clim([-0.5 0.5])
% subplot(212);bar(corrsmax)

% detect which components are not too correlated with dubious channels
badcomponents = abs(corrsmax) >= corr_threshold; 
find(badcomponents) % display the bad components
goodcomponents = abs(corrsmax) < corr_threshold;
sum(goodcomponents) % number of goodcomponents left thus far

% further detect bad ones in the goodcomponents by examing the weight
% proportion in A (mixing matrix)
proportion_threshold=0.6; % if any cahnnel weighted higher than 0.6 in A
chancomponents = zeros(size(icasig,1),1);
B = zeros(nchans,size(icasig,1)); % 208 channel x 208 components
for n = 1:size(icasig,1) % loop through all 208 components
    B(:,n)=A(:,n).^2 / sum(A(:,n).^2);
    chancomponents(n)=max(B(:,n));
    if chancomponents(n) > proportion_threshold
        display(['detect bad component - ' num2str(n)])
        goodcomponents(n)=0; % remove it from good components
    end
end
sum(goodcomponents)

% figure;
% subplot(211);imagesc(B);colormap('jet');clim([0 1]);ylabel('nchans');xlabel('component weights');
% subplot(212);bar(chancomponents)


% recombine data without bad components.
mixedsig=A(:,goodcomponents)*icasig(goodcomponents,:);

% % Examine before and after ica
% ScreenSize=get(0,'MonitorPositions');
% FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
% figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/2]);
% subplot(211);
% eeg_example=good_filtered_data;
% time_example=1:size(eeg_example,1);
% plot(time_example,eeg_example); 
% title('good-filtered-data');
% subplot(212);
% eeg_example=mixedsig';
% time_example=1:size(eeg_example,1);
% plot(time_example,eeg_example); 
% title('mixedsig');
% % selet two point on x axis to zoom in
% [x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% % get the proximate index
% string(x)
% ind1=round(x(1))
% ind2=round(x(2))
% subplot(211);
% eeg_example=good_filtered_data;
% time_example=1:size(eeg_example,1);
% plot(time_example(ind1:ind2),eeg_example(ind1:ind2,:)); 
% title('good-filtered-data');
% ylim([-50 50])
% subplot(212);
% eeg_example=mixedsig';
% time_example=1:size(eeg_example,1);
% plot(time_example(ind1:ind2),eeg_example(ind1:ind2,:)); 
% title('mixedsig');
% ylim([-50 50])

%% inserting the periphereal channels back in place with zeros
% need to use full 256 channel for mne template
preprocessed_eeg = zeros(256, size(mixedsig, 2));
i = 1;
for c = 1:256
    if ~ismember(c, ch_peripheral)
        preprocessed_eeg(c, :) = mixedsig(i, :);
        i = i + 1;
    end
end

clearvars A B corrsmax du icasig n W mixedsig corr_threshold good_filtered_data

%% get the original index of the dubious chans
open goodchanlocs;

ch_dubious=nan(1,length(dubious_chans));
for i = 1:length(dubious_chans)
    ch_dubious(i) = central_chanlocs(dubious_chans(i)).urchan;
end
dubious_chans; % 67    73    81   192   193   199   200   201   202   203   204   205   206   207   208
ch_dubious;    % 67    73    81   218   219   225   226   227   230   231   248   249   252   253   254

%% save the cleaned data
% filtered and cleaned data saved here
cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_cleaned/
save([num2str(subject_ID) '.mat'],'preprocessed_eeg','Fs','ch_dubious','ch_peripheral','ch_labels','chanlocs','subject_ID')
% Did not save 'goodchanlocs' (208 channels)

%% source localization and pca
% refer to 
% open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/loop_source_data.m
% load forward matrix
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/HEAD_model_ESCH
load('ESCH_ico_4_scale_0.05_depth_0.8.mat')
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
    vox = floor(source_fsaverage(i,:)); % change from ceil to floor,now we have 2 subcortical not mapped
    inds              = sub2ind([size(parcels)], vox(1), vox(2), vox(3));
    label             = parcels(inds); 
    source_labels(i) = label;
end

% find which label not mapped
roiNames_250(setdiff(1:463,unique(source_labels))) % 227 457 459
setdiff(1:463,unique(source_labels))

% inverse model
addpath ../../AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat, stat, reconstructed] = inversemodel(leadfield,'prctile',1);

% % check if all source_rr and leadfield for every subject is the same?
% cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/HEAD_model_ESCH
% load('ESCH_ico_4_scale_0.05_depth_0.8.mat')
% source_rr % 5124x3
% leadfield % 256x5124
% cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_spacing_DAON
% load('DAON_ico_4_scale_0.05_depth_0.8.mat')
% /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_spacing_JOHR
% load('JOHR_ico_4_scale_0.05_depth_0.8.mat')
% % confirmed that they are the same


source_data=inversemat*preprocessed_eeg;
EEG_recon=leadfield*source_data;

fra_eigenvalues=zeros(1,max(unique(source_labels)));
agr_source_data=[];
ave_source_coor=[];
ave_source_label=[];
tic
for sr=1:max(unique(source_labels))
    I=find(source_labels==sr);
    if ~isempty(I)
        [COEFF, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
        fra_eigenvalues(sr)=LATENT(1)/sum(LATENT);
        agr_source_data=[agr_source_data SCORE(:,1)];
        ave_source_coor=[ave_source_coor; mean(source_fsaverage(I,:),1)];
        ave_source_label=[ave_source_label; sr];
    end
end
toc % 66s

cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/organize_62stroke
% save([num2str(subject_ID) '.mat'],'fra_eigenvalues','agr_source_data','ave_source_coor','ave_source_label');

% remove the non-zeros fraction of eigenvalues
marked_fra_eigenvalues=fra_eigenvalues(ave_source_label);
bar(fra_eigenvalues)
bar(marked_fra_eigenvalues)

% create a boolean of subcortical rois
bool_subcorti=zeros(1,length(ave_source_label));
for i=1:length(ave_source_label)
    clear tmp
    tmp=ave_source_label(i);
    if ismember(tmp, scale250_subcortROIs)
        bool_subcorti(i)=1;
    end
end
sum(bool_subcorti)
ind=find(bool_subcorti);
ind
% remove subcortical roi
agr_source_data(:,ind)=[];
marked_fra_eigenvalues(ind)=[];
ave_source_coor(ind,:)=[];
ave_source_label(ind)=[];

corti_fra_eigenvalues=marked_fra_eigenvalues;
corti_ave_source_coor=ave_source_coor;
corti_ave_source_labl=ave_source_label;

cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source
save('corti_ave_source_coor.mat','corti_ave_source_coor');
save('corti_ave_source_labl.mat','corti_ave_source_labl');

save([num2str(subject_ID) '.mat'],'corti_fra_eigenvalues','agr_source_data','corti_ave_source_coor','corti_ave_source_labl');


%% hilbert
% refer to 
% open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/hilbert2cov.m
sampl_rate=Fs;
srnew = 200;
downsample = sampl_rate/srnew;
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

downsample_data=resample(double(agr_source_data),1,downsample,'Dimension',1);


open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/hilbert2cov.m
% filter
hilbert_dataCov=cell(1,5);
tic
for freq=1:5
    filterd_data = filter(filt_ds{freq},downsample_data);
    hilbertdata = hilbert(filterd_data');
    sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));
    sourceDataReal = sourceDataReal*(1/mean(abs(sourceDataReal(:)))); % normalize data
    hilbert_dataCov{freq} = cov(sourceDataReal');
end
toc

cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_source_hilbert
save([num2str(subject_ID) '.mat'],'hilbert_dataCov');

%% source coherence
% convert to complex covariance
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
Complex_dataCov_all=nan(5,448,448);
for freq=1:5
    Complex_dataCov_all(freq,:,:)=r2c(hilbert_dataCov{freq});
end
coh_all=nan(5,448,448);
for freq=1:5
    coh_all(freq,:,:)=normalizeCSD(squeeze(Complex_dataCov_all(freq,:,:)));
end

%% partial coherence AGL
% AGL variables
allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
n_Lambdas=length(allLambdas); % number of lambda values
min_LamdaIn = min(allLambdas);
% cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source
% save('corti_ave_source_coor.mat','corti_ave_source_coor');
% save('corti_ave_source_labl.mat','corti_ave_source_labl');
corti_ave_source_labl;
corti_ave_source_coor;
source_labels=corti_ave_source_labl; % every trial is the same, load one trial
source_coor=corti_ave_source_coor;
% penaltyselection for 1 subject
addpath(genpath('/../../AdaptiveGraphicalLassoforParCoh/Simulations/util'))
addpath(genpath('../../AdaptiveGraphicalLassoforParCoh/AGL'))
load('../../Virtual-Tractography/ForZhibin/processed_data/scale250_Connectome.mat');
SC=logical(fc(source_labels,source_labels));
sum(triu(SC,1),"all") % 4990 edges

hilbert_dataCov=cell(1,5);
datapermuted=cell(1,5);
datapermuted_cov=cell(1,5);
tic
for freq=1:5
    filterd_data = filter(filt_ds{freq},downsample_data);
    hilbertdata = hilbert(filterd_data');
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
    hilbert_dataCov{freq} = cov(sourceDataReal);
    for n=1:n_split
        datapermuted_cov{freq}(n,:,:) = cov([squeeze(datapermuted{freq}(n,:,:))]');
    end
end
toc % 9 s

cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_source_hilbert
save([num2str(subject_ID) '.mat'],'hilbert_dataCov','datapermuted_cov');

clear agr_source_data datareshaped downsample_data EEG_recon filterd_data filtered_data
clear hilbertdata preprocessed_eeg sourceDataReal

% penalty selection and fit precision
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
penalizationIn_op=nan(1,5);
penalizationOut_op=nan(1,5);
minDev_op3=nan(1,5);
X=nan(5,896,896)
for freq=1:5
    tic
    dataCovs_op=squeeze(datapermuted_cov{freq});
    [penalizationIn_op(freq),penalizationOut_op(freq),minDev_op3(freq)]=penaltyselection( ...
        SC,allLambdas,allLambdasOut,dataCovs_op);
    dataCov=hilbert_dataCov{freq};
    [X(freq,:,:)] = fitprecision( ...
        SC,penalizationIn_op(freq),penalizationOut_op(freq),min_LamdaIn,dataCov);
    toc 
    % 4 ensames each freq takes 964 s = 82 hours for 62 files
    % 2 ensames each freq takes 440 s = 37 hours for 62 files
end

