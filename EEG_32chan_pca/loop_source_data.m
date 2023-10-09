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
    vox = ceil(source_fsaverage(i,:));
    inds              = sub2ind([size(parcels)], vox(1), vox(2), vox(3));
    label             = parcels(inds); 
    source_labels(i) = label;
end

% load forward matrix
load('leadfield.mat');

% inverse model
addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat, stat, reconstructed] = inversemodel(leadfield,'prctile',1);
% leakage=inversemat*leadfield;

%% loop through all data
data_path= '../../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 
tic
corr_32_all=zeros(numSes,2,12,32);
for r=1:numSes
    clear runid
    runid = num2str(seeds(r,:));
    % data=load([data_path  'clean_' runid '.mat']);
    data=load([data_path  'clean_' runid '.mat'],'dataL','dataR');
    for subj=1:2
        if subj==1
            preprocessed_eeg=data.dataL;
        else
            preprocessed_eeg=data.dataR;
        end
        for tr=1:12
            EEG_ori=preprocessed_eeg{tr}(:,1:32)';
            source_data=inversemat*EEG_ori;
            EEG_recon=leadfield*source_data;

            for ch=1:32
                R=corrcoef(EEG_recon(ch,:)',EEG_ori(ch,:)');
                corr_32_all(r,subj,tr,ch)=R(1,2);
            end
        end
    end
end
toc
% 440 s 

for r=1:numSes
    clear runid data preprocessed_eeg EEG_ori source_data EEG_recon leakage
    runid = num2str(seeds(r,:));
    % data=load([data_path  'clean_' runid '.mat']);
    data=load([data_path  'clean_' runid '.mat'],'dataL','dataR');
    for subj=1:2
        if subj==1
            preprocessed_eeg=data.dataL
        else
            preprocessed_eeg=data.dataR
        end
        for tr=1:12
            clear EEG_ori source_data EEG_recon leakage
            tic
            EEG_ori=preprocessed_eeg{tr}(:,1:32)';
            source_data=inversemat*EEG_ori;
            EEG_recon=leadfield*source_data;

            fra_eigenvalues=zeros(1,max(unique(source_labels)));
            agr_source_data=[];
            ave_source_coor=[];
            ave_source_label=[];
            parfor sr=1:max(unique(source_labels))
                I=find(source_labels==sr);
                if ~isempty(I)
                    [COEFF, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
                    fra_eigenvalues(sr)=LATENT(1)/sum(LATENT);
                    if fra_eigenvalues(sr) > 0.5
                        agr_source_data=[agr_source_data SCORE(:,1)];
                        ave_source_coor=[ave_source_coor; mean(source_fsaverage(I,:),1)];
                        ave_source_label=[ave_source_label; sr];
                        save([data_path num2str(runid) '/subj' num2str(subj) '_tr_' num2str(tr)  '.mat'], ...
                            'fra_eigenvalues','agr_source_data','ave_source_coor','ave_source_label', ...
                            'leadkage','source_data');
                    end
                end
            end
            toc
        end
    end
end







