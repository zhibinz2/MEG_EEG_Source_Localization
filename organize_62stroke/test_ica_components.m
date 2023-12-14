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


%%  
for f=48
    tic
    cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_reorganized
    load([num2str(subj_files(f)) '.mat']);
    subject_ID=subj_files(f);
    display(['start processing subject file: ' num2str(subj_files(f)) '.mat']);

    % down-sample and z-score first


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

    % detect which components are not too correlated with any dubious channels
    badcomponents = abs(corrsmax) >= corr_threshold; 
    display(['bad components: ' num2str(find(badcomponents))] )% display the bad components
    goodcomponents = abs(corrsmax) < corr_threshold;
    display(['number of good components thus far: ' num2str(sum(goodcomponents))])  % number of goodcomponents left thus far
    ica_good_1(f)=sum(goodcomponents);

    % further detect bad ones in the goodcomponents by examing the weight
    % proportion in A (mixing matrix)
    proportion_threshold=0.6; % if any channel weighted higher than 0.6 in A of a component
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
    ica_good_2(f)=sum(goodcomponents);

    clear icasig A W 
    
    display(['Complete one file: ' num2str(subj_files(f)) '.mat ********************'])
    toc
end

save('ica_good_component_left.mat','ica_good_1','ica_good_2');

%% examine eeg quality
% 44, 48, 60 have very low number of component left
% 43.mat, 47.mat, 59.mat
cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_cleaned
figure;
for f=1:61
    load([num2str(f-1) '.mat'])
    plotdata=preprocessed_eeg(150:160,10000:11000)';
    plottime=1:size(plotdata,1);
    clf;
    plotx(plottime,plotdata);
    ylim([-100 100])
    title(num2str(f-1));
    subtitle(ica_good_2(f))
    pause(0.25)
end

figure;
imagesc(B);ylabel('channel');xlabel('components');colorbar