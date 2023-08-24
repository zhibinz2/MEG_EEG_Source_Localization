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
clearvars Data rData detrend_Data detrend_pad pad_hp pad_lp

%% set up channel info for topoplots (skip)
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat');
topoplot([1:32],chaninfo,'nosedir','+X','style','map','electrodes','off');colormap('hot')

cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/test_scripts/EEG_chan
load('egixyz.mat')
load('egichaninfo.mat')
load('leadfield_nn_rr.mat')
plot(egichaninfo.X,egichaninfo.Y,'r.')
X=egichaninfo.X; Y=egichaninfo.Y;
for ch=1:256
    egiXYinfo(ch).labels=ch_labels{ch};
    egiXYinfo(ch).X=X(ch);
    egiXYinfo(ch).Y=Y(ch);
    egitopoinfo(ch).urchan=ch;
end
topoplot([1:256],egiXYinfo,'nosedir','+X','style','map','electrodes','off');colormap('hot')

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
labeled_ch=[18 37 36 224 59 183 69 202 87 153 96 170 94 190 116 150 31 21 101 126];

x=x/100; y=y/120; z=z/110;
[theta,rho,z] = cart2pol(x,y,z);
[azimuth,elevation,r] = cart2sph(x,y,z);

for ch=1:256
    egitopoinfo(ch).labels=ch_labels{ch};
    egitopoinfo(ch).type=[];
    egitopoinfo(ch).theta=round(rad2deg(theta(ch)));
    egitopoinfo(ch).radius=rho(ch);
    egitopoinfo(ch).X=x(ch);
    egitopoinfo(ch).Y=y(ch);
    egitopoinfo(ch).Z=z(ch);
    egitopoinfo(ch).sph_theta=[];% round(rad2deg(azimuth(ch)));
    egitopoinfo(ch).sph_phi=[];%elevation(ch);
    egitopoinfo(ch).sph_radius=[];%r(ch);
    egitopoinfo(ch).urchan=ch;
    egitopoinfo(ch).ref=[];
    egitopoinfo(ch).impedance=[];
    egitopoinfo(ch).median_impedance=[];
end

topoplot([1:256],egitopoinfo,'nosedir','+X','style','map','electrodes','on');colormap('hot')
%% extract good chanlocs
cd /home/zhibinz2/Documents/GitHub/archive/EEG_StrokePatients_n61
chanlocs = readlocs('GSN-HydroCel-257.sfp');
figure;
subplot(121)
topoplot([1:257],chanlocs,'nosedir','+X','style','map','electrodes','on');
colormap('jet');colorbar;clim([1 257])
% remove bad chan
ch_bad=[241 242 243 238 239 240 ...
    244 245 246 247 251 256 91 102 111 120 133 145 165 174 187 199 208 216 229 233 237 236 235 234 ...
    232 228 217 209 200 188 175 166 156 146 134 121 112 103 92 82 255 250 257];
goodchanlocs=chanlocs;
goodchanlocs(ch_bad)=[];
subplot(122)
topoplot([1:208],goodchanlocs,'nosedir','+X','style','map','electrodes','on');
colormap('jet');colorbar;clim([1 257])

%%



