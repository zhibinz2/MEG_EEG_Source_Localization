clear
cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_reorganized
subj_files=[0:1:60];

% for the removal of 49 peripheral channels and the Cz (257-49=208) 
ch_peripheral_cz=[241 242 243 238 239 240 ...
    244 245 246 247 251 256 91 102 111 120 133 145 165 174 187 199 208 216 229 233 237 236 235 234 ...
    232 228 217 209 200 188 175 166 156 146 134 121 112 103 92 82 255 250 ...
    257];
ch_peripheral_cz=sort(ch_peripheral_cz);
% remove Cz channel from ch_peripheral_cz 
ch_peripheral=ch_peripheral_cz;
ch_peripheral(end)=[]; 

% original sampling frequency
Fs=1000;

% high pass (no paddings needed)
Hp = makefilter(Fs,0.25,0.01,6,20,0); % for keeping readiness potential
% Low pass (take a few seconds)
Lp = makefilter(Fs,50,51,6,20,0);  

% identify bad chan and epochs by standard deviation
var_threshold = 2.5;  % normalized variance threshold to reject trials.
chan_threshold = 2.5;  % to identify % (smaller values are stricter)
corr_threshold = 0.4;  % threshold for identifying whether an ICA component contains eye movement. (smaller values are stricter)

% organize channel info
chanlocs = load('/home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_streamline_DAON/chanlocs.mat');
chanlocs = chanlocs.chanlocs;
chanlocs(257)=[]; % Channel 257 is Cz, which is reference (values are all zeros)
% we have to keep only 256 channel to meet MNE's template

% Add the original sequence index
for ch=1:256 
    chanlocs(ch).urchan=ch;
end

% Create channel labels
ch_labels = cell(256,1);
for c = 1:256
    ch_labels{c}=num2str(c);
end
clear c ch

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
% ch_labels{257}='Cz';

% remove 49 perpheral channels (then to check channel index for eyechans later)
central_chanlocs=chanlocs;
central_chanlocs(ch_peripheral)=[];

% use "urchan" to get an array of the orginal index in central_chanlocs
urchan=extractfield(central_chanlocs,'urchan');

% set the eye channels for ICA
eyechans = [248 252 253 67 249 254 73 230 226 225 219 231 227 218]; eyechans=sort(eyechans);
% the above eyechans  67    73   218   219   225   226   227   230   231   248   249   252   253   254
% corresponds to      67    73   192   193   199   200   201   202   203   204   205   206   207   208 
% in the updated chanel index of central_chanlocs (icachans below)

% identify the updated indices in central_chanlocs for eyechan index 
eyechans_ind=nan(1,length(eyechans));
for c=1:length(eyechans)
    eyechans_ind(c)=find(urchan==eyechans(c));
end


%% Loop through the raw data, filter, clean and run ica

for f=1:61
    tic
    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_reorganized
    load([num2str(subj_files(f)) '.mat']);
    subject_ID=subj_files(f);
    display(['start processing subject file: ' num2str(subj_files(f)) '.mat']);

    % average re-reference
    rData=Data-ones(size(Data,1),1)*mean(Data,1);
    % remove 49 peripheral channels and cz reference
    rData(ch_peripheral_cz,:)=[];
    rData=rData';
    % detrend the EEG data (no padding needed)
    detrend_data=detrend(rData,1);
    % add paddings
    padding=zeros(round(size(detrend_data,1)/10), size(detrend_data,2));
    detrend_pad=cat(1,padding,detrend_data,padding);
    % high pass (no paddings needed)
    pad_hp=filtfilthd(Hp,detrend_pad);
    % Low pass (take a few seconds)
    pad_lp=filtfilthd(Lp,pad_hp);
    % remove paddings
    filtered_data=pad_lp((size(padding,1)+1):(size(padding,1)+size(detrend_data,1)),:);

    clearvars Data rData detrend_data detrend_pad pad_hp pad_lp padding
    % keep output filtered_data

    % organized into 1s epochs
    nepochs=floor(length(filtered_data)/Fs);
    nchans=size(filtered_data,2);
    epochdata=zeros(nepochs,Fs,nchans);
    for e = 1: nepochs
        epochdata(e,:,:)=filtered_data((1+(e-1)*Fs):(Fs+(e-1)*Fs),:);
    end
    clearvars filtered_data e 
    % keep output epochdata

    % identify bad chan and epochs by standard deviation
    eegstd=squeeze(std(epochdata,[],1));
    chanstd=sum(eegstd,2);
    epochstd=sum(eegstd,1);
    
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

    clearvars epochdata epochstd eegstd cahnstd good_epochdata n
    % keep output good_filtered_data

    % run ICA
    [icasig, A, W] = fastica(good_filtered_data');  % a few minutes

    % combine badchan into dubious chans with to compute correlation with
    % eye channels (icachans)
    dubious_chans =unique([badchan eyechans_ind]);

    % compute correlation with dubious chans (a few seconds)
    corrs=zeros(length(dubious_chans),size(icasig,1));
    for du = 1:length(dubious_chans)
        corrs(du, :)=corr(icasig',good_filtered_data(:,dubious_chans(du)));
    end
    corrsmax = max(corrs);

    % detect which components are not too correlated with dubious channels
    badcomponents = abs(corrsmax) >= corr_threshold; 
    display(['bad components: ' num2str(find(badcomponents))] )% display the bad components
    goodcomponents = abs(corrsmax) < corr_threshold;
    display(['number of good components thus far: ' num2str(sum(goodcomponents))])  % number of goodcomponents left thus far

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
    display(['good components left: ' num2str(sum(goodcomponents))])
    
    % recombine data without bad components.
    mixedsig=A(:,goodcomponents)*icasig(goodcomponents,:);  

    % inserting the periphereal channels back in place with zeros
    % need to use full 256 channel for mne template
    preprocessed_eeg = zeros(256, size(mixedsig, 2));
    i = 1;
    for c = 1:256
        if ~ismember(c, ch_peripheral)
            preprocessed_eeg(c, :) = mixedsig(i, :);
            i = i + 1;
        end
    end
    clearvars icasig A W corrsmax du  n  mixedsig good_filtered_data
    % output preprocessed_eeg

    % get the original index of the dubious chans
    ch_dubious=nan(1,length(dubious_chans));
    for i = 1:length(dubious_chans)
        ch_dubious(i) = central_chanlocs(dubious_chans(i)).urchan;
    end

    % save the cleaned data
    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_cleaned/
    save([num2str(subject_ID) '.mat'],'preprocessed_eeg','Fs','ch_dubious','ch_peripheral','ch_labels','chanlocs','subject_ID')
    
    display(['Complete one file: ' num2str(subj_files(f)) '.mat ********************'])
    toc
end

%% % source localization and pca to get cortical source data
% load forward matrix and source
clear
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/HEAD_model_ESCH
load('ESCH_ico_4_scale_0.05_depth_0.8.mat','leadfield','source_rr');
load('Lausanne2008_fsaverageDSsurf_60_125_250.mat','Brain', ...
    'roiNames_250','scale250_subcortROIs')
load('parcels.mat') % This is the labels
% label the sources
Vertex=Brain.Vertex;
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

% get coordinate and label for each aggregated ROIs
ave_source_coor=[];
ave_source_label=[];
for sr=1:max(unique(source_labels))
    I=find(source_labels==sr);
    if ~isempty(I)
        ave_source_coor=[ave_source_coor; mean(source_fsaverage(I,:),1)];
        ave_source_label=[ave_source_label; sr];
    end
end 

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
ind % use these subcortical indices to remove subcortical aggregated pca data
% remove subcortical roi
ave_source_coor(ind,:)=[];
ave_source_label(ind)=[];
corti_ave_source_coor=ave_source_coor;
corti_ave_source_labl=ave_source_label;
cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source
save('corti_ave_source_coor.mat','corti_ave_source_coor');
save('corti_ave_source_labl.mat','corti_ave_source_labl');

% get inverse model
addpath ../../AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat, stat, reconstructed] = inversemodel(leadfield,'prctile',1);



% source loclization and pca
subj_files=[0:1:60];
for f=1:61
    tic
    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_cleaned/
    load([num2str(subj_files(f)) '.mat'],'preprocessed_eeg','subject_ID', ...
        'chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');
    display(['start processing subject file: ' num2str(subj_files(f)) '.mat']);

    source_data=inversemat*preprocessed_eeg;

    fra_eigenvalues=zeros(1,max(unique(source_labels)));
    corti_source_data=[];
    for sr=1:max(unique(source_labels))
        I=find(source_labels==sr);
        if ~isempty(I)
            [COEFF, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
            fra_eigenvalues(sr)=LATENT(1)/sum(LATENT);
            corti_source_data=[corti_source_data SCORE(:,1)];
        end
    end % 66s

    % remove the non-zeros fraction of eigenvalues and subcortical rois
    corti_fra_eigenvalues=fra_eigenvalues(ave_source_label);
    
    % remove subcortical roi
    corti_source_data(:,ind)=[];

    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source
    save([num2str(subject_ID) '.mat'],'corti_fra_eigenvalues', ...
        'corti_source_data','corti_ave_source_coor','corti_ave_source_labl', ...
        'subject_ID','chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');

    display(['Complete one file: ' num2str(subj_files(f)) '.mat ********************'])
    toc

end

%% downsample filter hilbert -> coh and partial coh
% refer to 
% open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/hilbert2cov.m
clear
Fs=1000;
srnew = 200;
downsample = Fs/srnew;
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
load('../../Virtual-Tractography/ForZhibin/processed_data/scale250_Connectome.mat','fc');
SC=logical(fc(corti_ave_source_labl,corti_ave_source_labl));
sum(triu(SC,1),"all") % 4990 edges

% penalty selection and fit precision
subj_files=[0:1:60];
for f=3:15
    tic
    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source
    load([num2str(subj_files(f)) '.mat'],'corti_source_data','subject_ID', ...
        'chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');
    display(['start processing subject file: ' num2str(subj_files(f)) '.mat']);

    downsample_data=resample(double(corti_source_data),1,downsample,'Dimension',1);
    
    datapermuted=cell(1,5);
    stroke_Cov=cell(1,5); % coh
    stroke_Pcov=cell(1,5); % partial coh
    penalizationIn_op=nan(1,5);
    penalizationOut_op=nan(1,5);
    minDev_op=nan(1,5);
    for freq=1:5
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

    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_source_hilbert
    save([num2str(subject_ID) '.mat'],'stroke_Cov', 'stroke_Pcov', ...
        'penalizationIn_op','penalizationOut_op','minDev_op', ...
        'subject_ID','chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');

    display(['Complete one file: ' num2str(subj_files(f)) '.mat ********************'])
    toc
end