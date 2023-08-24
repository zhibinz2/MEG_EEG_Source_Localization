% Author : Alexandre Chalard, Ph.D
% Date : 02/10/22
% This script create epoch for EEG data for both the task and the rest

%% Setup
clear
Folder = 'D:\UCLA\bmyoung\SM_EeEeGee\Data\Preprocessed'; % Folder with preprocess .set data
FolderOut = 'D:\UCLA\bmyoung\SM_EeEeGee\Data\Epoch'; % Folder where to save the data
SubToProcess = {'DEVT'}; % Enter ID of the subject to process or "All" for all the subjects in the RawFolder
cd(FolderOut)

% 'BATR','BEAD','CHAA','HANS','HONL','KATG','LIUC','MEJC','OPRO','PATD','SAKS','SHIK','YOUB'

%% Epoch for Task data

%------------------------
% Files selections
%------------------------
ListF  = ListFSubDir(Folder,'EEG',{'REST'});
if ~strcmp(SubToProcess,'All')
    ListF = ListF(contains(ListF,SubToProcess));
end

[~,files] = cellfun(@fileparts,ListF,'UniformOutput',false);
spfiles   = squeeze(split(files,'_'));
Sub       = unique(spfiles(:,1));


for sub = 1 : numel(Sub)
    %--------------------------------
    % Loading and formatting data ---
    % -------------------------------
    
    uListF     = ListFSubDir(Folder,strcat(Sub{sub}),{'REST'});
    SpeedF     = ListFSubDir(fullfile('D:\UCLA\bmyoung\SM_EeEeGee\Data\Raw',Sub{sub}),'InitialSpeed');
    load(SpeedF{1})
    EpochTime = Speed{end}(1)+0.5; % Epoch 0.5s after movement
    DataSet = pop_loadset('filename',uListF);
    
    
    %----------------------------------------
    % Merging file in one & create epoch  ---
    % --------------------------------------
    
    
    DataSet = pop_epoch(DataSet,{},[-2 EpochTime]);
    ALLEEG = [];
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, DataSet, 0,'study',0);
    EEG = eeg_checkset(EEG);
    Merged = pop_mergeset( ALLEEG,1:length(EEG), 0);
    clear EEG
    EEG = Merged;
    
    
    
    %--------------------------------------
    % Automatic rejection of bad epochs ---
    % -----------------------------------
    
    % For automatic rejection see https://doi.org/10.1016/j.neuroimage.2019.06.046
    % Compute metrics described in the paper and reject the 15 % baddest
    % epoch
    
    for n = 1 : size(EEG.data,3)
        ep = EEG.data(:,:,n);
        % This parameter have been set up according to the default and
        % methods described in paper, feel free to play with it to see
        % the change
        settings.overallThresh = 30;
        settings.timeThresh = 15;
        settings.chanThresh = 15;
        [OHA(n),THV(n),CHV(n)] = GetEpochMetrics(ep,settings);
        
    end
    
    % Reject 10% of the worst epoch feel free to change it
    RejThreshold = 10;
    iOHA = find(OHA > prctile(OHA,100-RejThreshold)==1);
    iTHV = find(THV > prctile(THV,100-RejThreshold)==1);
    iCHV = find(CHV > prctile(CHV,100-RejThreshold)==1);
    Idx = unique([iOHA iTHV iCHV]);
    
    EEG = pop_rejepoch(EEG,Idx,0);
    
    data = eeglab2fieldtrip(EEG,'raw','none');
    
    cfg = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-2 -1.75];
    data                = ft_preprocessing(cfg,data);
    
    
    save(strcat(Sub{sub},'_Task.mat'),'data');
    
    clearvars -except ListF files Sub LabFreq Folder FolderOut SubToProcess sub f
    
end
clearvars -except Folder FolderOut SubToProcess


%% Epoch for REST for a non band pass REST data

%------------------------
% Files selections
%------------------------

ListF  = ListFSubDir('D:\UCLA\bmyoung\SM_EeEeGee\Data\Preprocessed','REST');
if ~strcmp(SubToProcess,'All')
    ListF = ListF(contains(ListF,SubToProcess));
end

for file = 1 : numel(ListF)
    
    EEG = pop_loadset('filename',ListF{file});
    EEG = eeg_regepochs(EEG);
    
    %--------------------------------------
    % Automatic rejection of bad epochs ---
    % -----------------------------------
    
    % For automatic rejection see https://doi.org/10.1016/j.neuroimage.2019.06.046
    % Compute metrics described in the paper and reject the 15 % baddest
    % epoch
    
    for n = 1 : size(EEG.data,3)
        ep = EEG.data(:,:,n);
        settings.overallThresh = 30;
        settings.timeThresh = 15;
        settings.chanThresh = 15;
        [OHA(n),THV(n),CHV(n)] = GetEpochMetrics(ep,settings);
        
    end
    
    % Reject 10% of the worst epoch feel free to change it
    RejThreshold = 10;
    iOHA = find(OHA > prctile(OHA,100-RejThreshold)==1);
    iTHV = find(THV > prctile(THV,100-RejThreshold)==1);
    iCHV = find(CHV > prctile(CHV,100-RejThreshold)==1);
    Idx = unique([iOHA iTHV iCHV]);
    
    EEG = pop_rejepoch(EEG,Idx,0);
    
    data = eeglab2fieldtrip(EEG,'raw','none');
    
    cfg = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [0 1];
    data                = ft_preprocessing(cfg,data);
    
    
    [~,fName] = fileparts(ListF{file});
    spName = split(fName,'_');
    
    
    save(strcat(spName{1},'_REST.mat'),'data');
    
    clearvars -except ListF files Sub LabFreq Folder FolderOut SubToProcess sub f
    
    
end
