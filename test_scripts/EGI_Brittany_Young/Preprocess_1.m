% Author : Alexandre Chalard, Ph.D
% Date : 02/08/22
% EEG preprocessing scripts for the rotor pursuit task

%% Setup
clear
RawFolder = 'D:\UCLA\bmyoung\SM_EeEeGee\Data\Raw'; % Folder with raw .mat data
FunFolder = 'D:\UCLA\bmyoung\SM_EeEeGee\RotorPursuit_Scripts\fun'; % Folder where the function are
FolderOut = 'D:\UCLA\bmyoung\SM_EeEeGee\Data\Preprocessed'; % Folder where to save the data
SubToProcess = {'DEVT'}; % Enter ID of the subject to process or "All" for all the subjects in the RawFolder
SubjectFile  = 'D:\UCLA\bmyoung\SM_EeEeGee\SubjectData.xlsx';

% Double Two variables with same name, choose the good one after loading
%DEVT_EEG_Block_1.mat

% 'HANS','HONL','KATG','LIUC','MEJC','OPRO','PATD','SAKS','SHIK','YOUB'

%% Making magic

%------------------------
% Files selections
%------------------------

ListF  = ListFSubDir(RawFolder,'EEG');
if ~strcmp(SubToProcess,'All')
    ListF = ListF(contains(ListF,SubToProcess));
end

addpath(FunFolder)
eeglab nogui
RemoveChan = [67 73 82 91 92 102 103 111 112 120 121 133 134 145 146 156 165 ...
    166 174 175 187 188 199 200 208 209 216 217 218 219 225 226 227 228 229 230 ...
    231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 ...
    250 251 252 253 254 255 256 257];
cd(FolderOut)

tab = readtable(SubjectFile);


for f = 1 : numel(ListF)
    
    %--------------------------------
    % Loading and formatting data ---
    % -------------------------------
    
    display(strcat('Loading file : ',ListF{f}));
    load(ListF{f})
    [~,origname] = fileparts(ListF{f});
    spName   = split(origname,'_');
    name     = spName{1};
    vars = who();
    TF = contains(vars, name);
    NameVar  = vars(TF);
    
    if length(NameVar) > 1
        error("Two variables with same name, choose the good one")
        return
    end


% IF multiple EEG recordings are within the raw .mat that is loaded, select
% the correct file within NameVar (in the case of DEVT it is the file with
% index 2) and then select the remainder of the script and run the rest of
% the code. 

% If there are additional files from the affected subject to run afer the
% affected file, change the indices in the "for f = 1 : numel(ListF)" line
% (e.g. "for f=2:numel(ListF)") before re-running the code with that
% subject ID specifically specified. -- Don't forget to change the indices
% back to 1:numel(ListF) after running. 

% Remember to change the NameVar index back to 1 in the code below before
% running the script on future subjects. 

    Data = evalin('base',NameVar{1}); % Change the number inside NameVar to match the good recording

    Events = cell2mat(evt_255_DINs(4,ismember(evt_255_DINs(1,:),'D254')));
    chanlocs = readlocs('D:\UCLA\bmyoung\SM_EeEeGee\RotorPursuit_Scripts\fun\GSN-HydroCel-257.sfp');
    
    for ch = 1 : length(chanlocs)
        chanlocs(ch).type = 'EEG';
    end
    
    subID = spName{1};

    %------------------------------------------------
    %-- If Rest cut from white square to white square
    %-------------------------------------------------
    
    if any(contains(spName,'REST'))
        Data = Data(:,Events(1):Events(end));
    end
    
    
    
    samplingRate = EEGSamplingRate;
    EEG = pop_importdata('setname','DoubleArrow','data','Data','dataformat','array',...
        'nbchan',257,'chanlocs',chanlocs,'srate',samplingRate);
    
    
    %--------------------------------
    %-- Flip data if need (Left side)
    %--------------------------------
    
    side = table2cell(tab(ismember(tab.Study_ID,subID),8));
    if isempty(side)
        warning('Side not found check excel file');
        return
    end
    % If left affected or task performed on left side flip EEG data
    if strcmp(side,'L')
        EEG.data = Flip_RL(EEG.data);
    end
    

    %--------------------------------
    % Processing events for Task  ---
    % -------------------------------
    
    if ~any(contains(spName,'REST'))
        WacomFile = fullfile(RawFolder,spName{1},strcat(spName{1},'_Wacom_Block_',spName{4},'.mat')); % Load Data for event detection
        Cond = table2cell(tab(ismember(tab.Study_ID,subID),3));
        if strcmp(Cond,'Task')
            EventToRemove = findBadTrials(WacomFile); % Find if very bad trials Perf < 5 %
        else % If movement attempt do not remove data
            EventToRemove = [];
        end
        event = cell(length(Events),2);
        for i = 1 : length(event)
            event{i,1} = 'Rotor';
            event{i,2} = Events(i)./1000; % Convert second
        end
        if ~isempty(EventToRemove)
            event(EventToRemove,:) = []; % Remove if bad trials
        end
        EEG = pop_importevent(EEG,'event',event,'fields',{'type' 'latency'});
        EEG.eventRemoved = EventToRemove;
    end
    
    %-----------------------------------------------------
    % Filtering + Elliptical Notch + Resampling 200Hz ----
    % ----------------------------------------------------
    
    EEG = pop_resample(EEG,200);
    
    display('60Hz Notch filtering..')
    match    = 'stopband';  % Band to match exactly
    h        = fdesign.bandstop(57, 57.9, 61.1, 62,1, 100, 1, 200);
    Hd_notch = design(h, 'butter', 'MatchExactly', match);
    tmpData  = filtfilthd(Hd_notch,EEG.data');
    EEG.data  = tmpData';
    clear tmpData
    
    EEG = pop_select(EEG,'nochannel',RemoveChan);
    EEG = pop_reref(EEG,[]);
    
    % Remove baseline
    EEG.data = rmbase(EEG.data);
    
    originalEEG = EEG;
    
    % Clean_rawdata, highpass filter, line noise removal, and bad channels
    % detection
    EEG = pop_clean_rawdata(EEG,'FlatlineCriterion',5,'ChannelCriterion',0.8,...
        'LineNoiseCriterion',4,'Highpass',[0.5 1] ,...
        'BurstCriterion','off','WindowCriterion','off','BurstRejection','off',...
        'Distance','Euclidian','WindowCriterionTolerances','off' );
    
    % Interpolate bad channels
    EEG = pop_interp(EEG,originalEEG.chanlocs,'spherical');
    
    
    %--------------------------------
    % Automatic ICA rejection    ---
    % -------------------------------
    
    % ICA decomposition
    EEG = pop_runica(EEG,'icatype','runica','verbose','off');
    try
        % ICA rejection using ICLabel
        EEG = pop_iclabel(EEG,'default');
        EEG = pop_icflag(EEG,[0 0; 0.8 1; 0.8 1; 0.8 1; 0 0; 0 0; 0 0]); % Remove  Muscle eye component
        EEG = pop_subcomp(EEG, []);
    catch
        % Remove with ICAlabel not running')
    end
    
    
    SaveName = strcat(origname);
    pop_saveset(EEG,'filename',SaveName);
    
    clearvars -except RawFolder FunFolder ListF RemoveChan FolderOut SubToProcess SubjectFile tab
end



