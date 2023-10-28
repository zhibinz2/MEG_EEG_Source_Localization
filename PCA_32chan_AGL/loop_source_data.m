%% source allignment and load some variables
% cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
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
    vox = floor(source_fsaverage(i,:)); % change from ceil to floor,now we have 2 subcortical not mapped
    inds              = sub2ind([size(parcels)], vox(1), vox(2), vox(3));
    label             = parcels(inds); 
    source_labels(i) = label;
end

% find which label not mapped
roiNames_250(setdiff(1:463,unique(source_labels))) % 225 227

% load forward matrix
load('leadfield.mat');

% inverse model
addpath ../../AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat, stat, reconstructed] = inversemodel(leadfield,'prctile',1);

data_path= '../../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 

cd ../../Cleaned_data


%% loop throught all source data to compute pca and aggregate (all ROIs)
% cd ../../Cleaned_data/source_data
% for ses=1:12
%     mkdir(num2str(seeds(ses,:)));
% end

cd ../../MEG_EEG_Source_Localization/PCA_32chan_AGL
for r=8%1:numSes
    vars = ({'runid', 'data preprocessed_eeg', 'EEG_ori', 'source_data', 'EEG_recon', 'leakage'});
    clear(vars{:});
    runid = num2str(seeds(r,:));
    data=load([data_path  'clean_' runid '.mat'],'dataL','dataR');
    for subj=1:2
        if subj==1
            preprocessed_eeg=data.dataL;
        else
            preprocessed_eeg=data.dataR;
        end
        for tr=1:12
            vars = ({'EEG_ori', 'source_data', 'EEG_recon', 'leakage'});
            clear(vars{:});
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
                    agr_source_data=[agr_source_data SCORE(:,1)];
                    ave_source_coor=[ave_source_coor; mean(source_fsaverage(I,:),1)];
                    ave_source_label=[ave_source_label; sr];
                end
            end
            save([data_path 'source_data/' num2str(runid) '/subj' num2str(subj) '_tr_' num2str(tr)  '.mat'], ...
                            'fra_eigenvalues','agr_source_data','ave_source_coor','ave_source_label');
            toc
        end
    end
end
% 10 h


%% Remove 16 subcortical and "zeros marked" sources not mapped and saved the aggreated source data
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/cortical_source_data
for ses=1:12
    mkdir(num2str(seeds(ses,:)));
end
load('/home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/Lausanne2008_fsaverageDSsurf_60_125_250.mat')
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
            % remove subcortical ROIs
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
% 40 min

