%% pick the raw eeg data file from a subject
% I picked Subject DAON, EEG code EPHD21 for testing, becasue we have this subject's 
% raw EEG data shared by Cramer, and the lesion mask, and the processed EEG
% data EPHD21
clear
% raw EEG file from Cramer
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_streamline/
cd ../../archive/EEG_StrokePatients_n61/
subject_ID='DAON'
load([subject_ID '_EPHD21_20151116_1200.mat'])
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_streamline/

%% average re-reference
rData=Data-mean(Data,2)*ones(1,size(Data,2));

%% remove bad channels
% remove 49 channels
ch_bad=[241 242 243 238 239 240 ...
    244 245 246 247 251 256 91 102 111 120 133 145 165 174 187 199 208 216 229 233 237 236 235 234 ...
    232 228 217 209 200 188 175 166 156 146 134 121 112 103 92 82 255 250 257];
rData(ch_bad,:)=[];
rData=rData';

%% detrend the EEG data (no padding needed)
detrend_data=detrend(rData,1);

%% add paddings
padding=zeros(round(size(detrend_data,1)/10), size(detrend_data,2));
detrend_pad=cat(1,padding,detrend_data,padding);

%% High pass
% addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab_zhibin/EEG/hnl'))
Fs=1000;
% high pass (no paddings needed)
Hd = makefilter(Fs,0.25,0.01,6,20,0); % for keeping readiness potential
pad_hp=filtfilthd(Hd,detrend_pad);

%% Low pass
Hd = makefilter(Fs,50,51,6,20,0);  
pad_lp=filtfilthd(Hd,pad_hp);

%% remove paddings
filtered_data=pad_lp((size(padding,1)+1):(size(padding,1)+size(detrend_data,1)),:);

%% clear out memory before proceeding to next section
clearvars Data rData detrend_data detrend_pad pad_hp pad_lp
clearvars padding Hd

%% organized into 1s epochs
nepochs=floor(length(filtered_data)/Fs);
nchans=size(filtered_data,2);
epochdata=zeros(nepochs,Fs,nchans);
for e = 1: nepochs
    epochdata(e,:,:)=filtered_data((1+(e-1)*Fs):(Fs+(e-1)*Fs),:);
end
clearvars e

%% identify bad chan and epochs by standard deviation (skip)
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
clearvars chan_threshold var_threshold

%% extract good chanlocs
chanlocs = load('chanlocs.mat');
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

% remove ch_bad to check channel index for eyechans
goodchanlocs=chanlocs;
for ch=1:256
    goodchanlocs(ch).labels=ch_labels{ch};
end
goodchanlocs(ch_bad)=[];

% remove last invalid channel from ch_bad
ch_bad=[241 242 243 238 239 240 ...
    244 245 246 247 251 256 91 102 111 120 133 145 165 174 187 199 208 216 229 233 237 236 235 234 ...
    232 228 217 209 200 188 175 166 156 146 134 121 112 103 92 82 255 250];
ch_bad=sort(ch_bad)
ch_bad

clearvars c ch chanloc

%% run ICA
% addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab_zhibin/EEG/hnl'))
[icasig, A, W] = fastica(good_filtered_data');

% Calculate Correlation with dubious chans
eyechans = [248 252 253 219 231 227 218 248 252 253 67 249 254 73]; 
icachans = [204 206 207 193 203 201 192 204 206 207 67 205 208 73];
% corresponding to chanel label 248 252 253 219 231 227 218 248 252 253 67 249 254 73
dubious_chans =unique([badchan icachans]);
corrs=zeros(length(dubious_chans),size(icasig,1));
for du = 1:length(dubious_chans)
    corrs(du, :)=corr(icasig',good_filtered_data(:,dubious_chans(du)));
end
corrsmax = max(corrs);

% detect which components are not too correlated with eye channels
goodcomponents = abs(corrsmax) < corr_threshold;
chancomponents = zeros(size(icasig,1),1);
B = zeros(nchans,size(icasig,1));
for n = 1:size(icasig,1)
    B(:,n)=A(:,n).^2 / sum(A(:,n).^2);
    chancomponents(n)=max(B(:,n));
    if chancomponents(n) > 0.7
        goodcomponents(n)=0;
    end
end

% recombine data without bad components.
mixedsig=A(:,goodcomponents)*icasig(goodcomponents,:);

%% inserting the bad channels back in place with zeros
preprocessed_eeg = zeros(256, size(mixedsig, 2));
i = 1;
for c = 1:256
    if ~ismember(c, ch_bad)
        preprocessed_eeg(c, :) = mixedsig(i, :);
        i = i + 1;
    end
end

clearvars A B corrsmax du icasig n W mixedsig corr_threshold good_filtered_data

%% need to map dubious chan index to the original 256 index
goodchanlocs;

ch_dubious=nan(1,length(dubious_chans))
for i = 1:length(dubious_chans)
    ch_dubious(i) = goodchanlocs(dubious_chans(i)).urchan;
end

ch_dubious;

%% save
save('preprocessed_eeg.mat','preprocessed_eeg','Fs','ch_dubious','ch_bad','ch_labels','chanlocs','subject_ID')
