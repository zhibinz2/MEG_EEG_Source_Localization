% I picked Subject DAON, EEG code EPHD21 for testing, becasue we have this subject's 
% raw EEG data shared by Cramer, and the lesion mask, and the processed EEG
% data EPHD21
clear

% processed EEG data
% cd ../../archive/STROKE/CRAMER_directory/subacute_STROKE/EPHD21/EEG/

% lesion mask
% cd ../../archive/STROKE/lesion_masks/30_DAON_final_strokemask_MNI_flipped_bin.nii.gz

% raw EEG file from Cramer
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/test_scripts/EEG_processing
cd ../../../archive/EEG_StrokePatients_n61/
load('DAON_EPHD21_20151116_1200.mat')


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

%% Plot and see
% ScreenSize=get(0,'MonitorPositions');
% FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
% figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/4]);
% eeg_example=detrend_data;
% time_example=1:size(eeg_example,1);
% plot(time_example,eeg_example); 
% ylim([-1000 1000])
% title('detrend data');

%% add paddings
padding=zeros(round(size(detrend_data,1)/10), size(detrend_data,2));
detrend_pad=cat(1,padding,detrend_data,padding);

%% Plot and see
% ScreenSize=get(0,'MonitorPositions');
% FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
% figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/4]);
% eeg_example=detrend_pad;
% time_example=1:size(eeg_example,1);
% plot(time_example,eeg_example); 
% title('detrend-pad');


%% High pass
% addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab_zhibin/EEG/hnl'))

Fs=1000;

% high pass (no paddings needed)
Hd = makefilter(Fs,0.25,0.01,6,20,0); % for keeping readiness potential
pad_hp=filtfilthd(Hd,detrend_pad);

%% Plot and see
% ScreenSize=get(0,'MonitorPositions');
% FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
% figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/4]);
% eeg_example=pad_hp;
% time_example=1:size(eeg_example,1);
% plot(time_example,eeg_example); 
% title('pad-hp');

%% Low pass
Hd = makefilter(Fs,50,51,6,20,0);  
pad_lp=filtfilthd(Hd,pad_hp);

%% Plot and see
% ScreenSize=get(0,'MonitorPositions');
% FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
% figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/4]);
% eeg_example=pad_lp;
% time_example=1:size(eeg_example,1);
% plot(time_example,eeg_example); 
% title('pad-lp');

%% remove paddings
filtered_data=pad_lp((size(padding,1)+1):(size(padding,1)+size(detrend_data,1)),:);

%% Plot and see
% ScreenSize=get(0,'MonitorPositions');
% FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
% figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/4]);
% eeg_example=filtered_data;
% time_example=1:size(eeg_example,1);
% plot(time_example,eeg_example); 
% title('filtered-data');

%% clear out memory before proceeding to next section
clearvars Data rData detrend_data detrend_pad pad_hp pad_lp

%% organized into 1s epochs
nepochs=floor(length(filtered_data)/Fs);
nchans=size(filtered_data,2);
epochdata=zeros(nepochs,Fs,nchans);
for e = 1: nepochs
    epochdata(e,:,:)=filtered_data((1+(e-1)*Fs):(Fs+(e-1)*Fs),:);
end
% % pick a epoch to examine
% e =100
% plot(1:Fs,squeeze(epochdata(e,:,:)))
%% identify bad chan and epochs by standard deviation (skip)
eegstd=squeeze(std(epochdata,[],1));
chanstd=sum(eegstd,2);
epochstd=sum(eegstd,1);

var_threshold = 2.5;  % normalized variance threshold to reject trials.
chan_threshold = 2.5;  % to identify 
% (smaller values are stricter)
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


% ScreenSize=get(0,'MonitorPositions');
% FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
% figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/4]);
% eeg_example=good_filtered_data;
% time_example=1:size(eeg_example,1);
% plot(time_example,eeg_example); 
% title('good-filtered-data');


%% set up channel info for topoplots (skip)
% cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
% load('chaninfo.mat');
% topoplot([1:32],chaninfo,'nosedir','+X','style','map','electrodes','off');colormap('hot')
% 
% cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/test_scripts/EEG_chan
% load('egixyz.mat')
% load('egichaninfo.mat')
% load('leadfield_nn_rr.mat')
% plot(egichaninfo.X,egichaninfo.Y,'r.')
% X=egichaninfo.X; Y=egichaninfo.Y;
% for ch=1:256
%     egiXYinfo(ch).labels=ch_labels{ch};
%     egiXYinfo(ch).X=X(ch);
%     egiXYinfo(ch).Y=Y(ch);
%     egitopoinfo(ch).urchan=ch;
% end
% topoplot([1:256],egiXYinfo,'nosedir','+X','style','map','electrodes','off');colormap('hot')
% 
% % Create channel labels
% ch_labels = cell(256,1);
% for c = 1:256
%     ch_labels{c}=num2str(c);
% end
% ch_labels{18}='Fp2'; ch_labels{37}='Fp1';
% ch_labels{36}='F3'; ch_labels{224}='F4';
% ch_labels{59}='C3'; ch_labels{183}='C4';
% ch_labels{69}='T7'; ch_labels{202}='T8';
% ch_labels{87}='P3'; ch_labels{153}='P4';
% ch_labels{96}='P7'; ch_labels{170}='P8';
% ch_labels{94}='LM'; ch_labels{190}='RM';
% ch_labels{116}='O1'; ch_labels{150}='O2';
% ch_labels{31}='NAS'; ch_labels{21}='Fz'; ch_labels{101}='Pz'; ch_labels{126}='Oz';
% labeled_ch=[18 37 36 224 59 183 69 202 87 153 96 170 94 190 116 150 31 21 101 126];
% 
% x=x/100; y=y/120; z=z/110;
% [theta,rho,z] = cart2pol(x,y,z);
% [azimuth,elevation,r] = cart2sph(x,y,z);
% 
% for ch=1:256
%     egitopoinfo(ch).labels=ch_labels{ch};
%     egitopoinfo(ch).type=[];
%     egitopoinfo(ch).theta=round(rad2deg(theta(ch)));
%     egitopoinfo(ch).radius=rho(ch);
%     egitopoinfo(ch).X=x(ch);
%     egitopoinfo(ch).Y=y(ch);
%     egitopoinfo(ch).Z=z(ch);
%     egitopoinfo(ch).sph_theta=[];% round(rad2deg(azimuth(ch)));
%     egitopoinfo(ch).sph_phi=[];%elevation(ch);
%     egitopoinfo(ch).sph_radius=[];%r(ch);
%     egitopoinfo(ch).urchan=ch;
%     egitopoinfo(ch).ref=[];
%     egitopoinfo(ch).impedance=[];
%     egitopoinfo(ch).median_impedance=[];
% end
% 
% topoplot([1:256],egitopoinfo,'nosedir','+X','style','map','electrodes','on');colormap('hot')
%% extract good chanlocs
% cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/test_scripts/EGI_Brittany_Young
cd /home/zhibinz2/Documents/GitHub/archive/EEG_StrokePatients_n61
chanlocs = readlocs('GSN-HydroCel-257.sfp');
% figure;
% subplot(121)
% topoplot([1:257],chanlocs,'nosedir','+X','style','map','electrodes','labels');
% colormap('jet');colorbar;clim([1 257])

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
ch_labels{31}='NAS'; ch_labels{21}='Fz'; ch_labels{101}='Pz'; ch_labels{126}='Oz';
% labeled_ch=[18 37 36 224 59 183 69 202 87 153 96 170 94 190 116 150 31 21 101 126];

goodchanlocs=chanlocs;
for ch=1:256
    goodchanlocs(ch).labels=ch_labels{ch};
end

% remove bad chan of mine
goodchanlocs(ch_bad)=[];
% subplot(122)
% topoplot([1:208],goodchanlocs,'nosedir','+X','style','map','electrodes','labels');
% colormap('jet');colorbar;clim([1 257])


%% run ICA
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/test_scripts/EEG_processing
% addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab_zhibin/EEG/hnl'))
[icasig, A, W] = fastica(good_filtered_data');

% % Plot all ICA component
% for i=1:size(A,2)
%     SqA(i)=sumsqr(A(:,i));
% end
% figure;
% plot(1:size(A,2),SqA,'ro');ylabel('sum of square of column in A');xlabel('ICs'); 
% title('all components');
% 
% [B,I]=sort(SqA,'descend');
% ComponentsExam=I(1:10);

% Calculate Correlation with dubious chans
eyechans = [202 10 18 37 46 204]; % chanel label 230 10 18 37 46 248
dubious_chans =unique([badchan eyechans]);
corrs=zeros(length(dubious_chans),size(icasig,1));
for du = 1:length(dubious_chans)
    corrs(du, :)=corr(icasig',good_filtered_data(:,dubious_chans(du)));
end
% imagesc(coors);colormap('jet');colorbar;xlabel('components');ylabel('dubious chans'); clim([-0.5 0.5])
corrsmax = max(corrs);
% bar(corrmax)

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
% bar(chancomponents)

% recombine data without bad components.
mixedsig=A(:,goodcomponents)*icasig(goodcomponents,:);

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

%% clear out memory before proceeding to next section
clearvars eeg_example detrend_data epochdata good_epochdata good_filtered_data

%% save
save('icadata.mat','mixedsig','Fs',"goodchan",'badchan')

%% save changinfo
% good_labels=ch_labels;
% good_labels(ch_bad(1:end-1))=[];
% for i = 1:length(goodchanlocs)
%     X(i)=goodchanlocs(i).X;
%     Y(i)=goodchanlocs(i).Y;
%     Z(i)=goodchanlocs(i).Z;