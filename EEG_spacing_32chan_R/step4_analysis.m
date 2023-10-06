%%
clear
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_spacing_32chan_R
subject_ID='20220804'; 
scales5=[0.01, 0.05, 0.1, 0.3, 0.5]; 
% methods={'MNE','dSPM','sLORETA','eLORETA'};
methods={'MNE','sLORETA','eLORETA'};
% load([subject_ID '_scale_' num2str(scale) '_depth_' num2str(depth) '.mat']);

%% (Color Scheme) 
red   = [1 0 0];
pink  = [1 0.65 0.75];
black = [0 0 0];
blue  = [0 0 1];
darkgreen = [0 0.5 0];
deepyellow  = [1 0.8 0.2];
matlab_blue=[0 0.4470 0.7410];
matlab_orange=[0.8500 0.3250 0.0980];
matlab_purple=[0.4940 0.1840 0.5560];
% combine colors
% colors9=[red;pink;black;blue;darkgreen;deepyellow;matlab_blue;matlab_orange;matlab_purple];
% colors4=[matlab_blue;matlab_orange;darkgreen];
colors5=[matlab_purple;matlab_blue;matlab_orange;darkgreen; deepyellow];
%% examine 32 chan
scales={'0.01', '0.05', '0.1', '0.3', '0.5'};
depths={'0.8'};
i=1;d=2;
icos={'4'}; % icos={'2','3','4','5'}
for m=1:length(methods)
    figure;
    legends=cell(1,length(icos));
    for s =1:length(scales)
        hold on;
        filename=[subject_ID '_method_' methods{m} '_ico_' icos{i} '_scale_' scales{s} '_depth_' depths{d} '.mat'];
        load(filename);
        hold on;
        plot(1:length(corrcoef_diag),corrcoef_diag,'.', 'markersize', 5, 'color',colors5(s,:));
        legends{s}=['scale ' scales{s}];
        clear corrcoef_diag
    end
    legend(legends,'location','southeast')
    xlabel('channel');ylabel('corrcoef');
    title(['Correlation between original and reconstructed EEG - method : ' methods{m}]);
    subtitle(['subject-ID: ' subject_ID '-ico-' icos{i} '-depth-' depths{d}]);
end
%% examine evarage
for m=1:length(methods)
    corrcoef_var=nan(1,length(scales));
    for s = 1:length(scales)
        legends=cell(1,length(scales));
        hold on;
        filename=[subject_ID '_method_' methods{m} '_ico_' icos{i} '_scale_' scales{s} '_depth_' depths{d} '.mat'];
        load(filename);
        corrcoef_var(s)=mean(corrcoef_diag);
        clear corrcoef_diag
    end
    
    figure;
    scatter(scales5,corrcoef_var);
    xlabel('scale');ylabel('corrcoef-ave');
    title(['Correlation between original and reconstructed EEG - method : ' methods{m}]);
    subtitle(['subject-ID: ' subject_ID '-ico-' icos{i} '-depth-' depths{d}]);
    % xlim([3.5 4.5]); ylim([0.8 0.88]);
    grid on
end

%% examine leakage by muliplying invmat and gain
leakage=invmat*leadfield;
figure;
imagesc(leakage)
colorbar()
colormap('jet')
clim([-0.005 0.005])
title(['leakage: invmat x gain  ------   Method: ' methods{1}])
subtitle(['subject-ID: ' subject_ID '-ico-' icos{1} '-depth-' depths{1}]);
title('0.05')

%% tryout pca
clear
load hald
[COEFF, SCORE, LATENT] = pca(ingredients,'Centered',false);
recon_ingredients= SCORE(:,1:3)*COEFF(:,1:3)';
%% tryout svd
A = [1 2; 3 4; 5 6; 7 8]
[U,S,V] = svd(A)
U*S(:,1)*V(:,1)'

%% load source data after beamformer
load('20220713_method_MNE_ico_4_scale_0.05_depth_0.8.mat')
% pickle = py.importlib.import_module('pickle');
% data=h5read('stc.fif-stc.h5','/home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_spacing_32chan')
tic
load('20220713_method_MNE_ico_4_scale_0.05_depth_0.8_st.data1.mat')
load('20220713_method_MNE_ico_4_scale_0.05_depth_0.8_st.data2.mat')
load('20220713_method_MNE_ico_4_scale_0.05_depth_0.8_st.data3.mat')
load('20220713_method_MNE_ico_4_scale_0.05_depth_0.8_st.data4.mat')
load('20220713_method_MNE_ico_4_scale_0.05_depth_0.8_st.data5.mat')
toc #1min

% matlab crashed excuting the following
% stc_data=cat(1,stc_data1,stc_data2,stc_data3,stc_data4,stc_data5)

%% try inversemat.m 
cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat, stat,reconstructed] = inversemodel(leadfield,'prctile',75);
leakage=inversemat*leadfield;
figure;
imagesc(leakage)
colorbar()
colormap('jet')
clim([-0.005 0.005])
title('75')
% try different percentage and examine the correlation again

%% normalize leakage by rows
clear
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_spacing_32chan
load('20220713_method_MNE_ico_4_scale_0.5_depth_0.8.mat')
leakage=invmat*leadfield;
figure;
subplot(121)
imagesc(leakage)
colorbar()
colormap('jet')
clim([-0.005 0.005])
title('invmat*leadfield')
subplot(122)
imagesc(leakage./(diag(leakage)*ones(1,size(leakage,2))));
colorbar()
colormap('jet')
clim([-1 1])
title('invmat*leadfield (normalized by diagonal value row by row)')
%% compare sim_cov with 'auto' 'shrunk'
scales={'0.01', '0.05', '0.1', '0.3', '0.5'};
i=1;d=2; m=1;
icos={'4'}; % icos={'2','3','4','5'}

legends=cell(1,2);
figure;
clf

filename=[subject_ID '_method_' methods{m} '_ico_4_scale_' scales{2} '_depth_0.8' '.mat'];
load(filename);
hold on;
plot(1:length(corrcoef_diag),corrcoef_diag,'b.', 'markersize', 5);
legends{1}=['Regularization scale 0.05 we picked on the whole data'];
clear corrcoef_diag

hold on;
% load('20220713_method_MNE_ico_4_scale_shrunk_depth_0.8.mat')
load('20220713_method_MNE_ico_4_scale_auto_depth_0.8.mat')
plot(1:length(corrcoef_diag),corrcoef_diag,'R.', 'markersize', 5);
legends{2}=['MNE''s built-in ''auto'' regularization on a random 200ms data'];

legend(legends,'location','southeast')
xlabel('channel');ylabel('corrcoef');
title(['Correlation between original and reconstructed EEG']);
subtitle(['subject-ID: ' subject_ID '-method-MNE-ico-4-depth-0.08']);

%% compare sim_cov with scale of 0.01 and 0.99

scales={'0.01', '0.99'};
i=1;d=2; m=1;
icos={'4'}; % icos={'2','3','4','5'}

legends=cell(1,2);
figure;
clf
load('20220713_method_MNE_ico_4_scale_0.01_depth_0.8.mat')
hold on;
plot(1:length(corrcoef_diag),corrcoef_diag,'b.', 'markersize', 5);
legends{1}=['Regularization scale 0.01 we picked on the whole data'];
clear corrcoef_diag

hold on;
load('20220713_method_MNE_ico_4_scale_0.99_depth_0.8.mat')
plot(1:length(corrcoef_diag),corrcoef_diag,'R.', 'markersize', 5);
legends{2}=['Regularization scale 0.99 we picked on the whole data'];

legend(legends,'location','southeast')
xlabel('channel');ylabel('corrcoef');
title(['Correlation between original and reconstructed EEG']);
subtitle(['subject-ID: ' subject_ID '-method-MNE-ico-4-depth-0.08']);

diag(sim_cov_mat)
imagesc(sim_cov_mat);colorbar;
%% runWeightedL2norm
cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util/%runWeightedL2norm.m
L=leadfield;
leakage=[squeeze(allOpInv(5,:,:))]'*leadfield;

%% test the performance of inversemat.m 
load('preprocessed_eeg.mat')
load('20220804_scale_0.05_depth_0.8.mat')

cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
prc = [1:99];
% prc = [1 5 10 25 50 75];
% corr_all32=nan(6,32);
corr_all=nan(1,length(prc));
leakage_all=nan(1,length(prc));
for i=1:length(prc)
    tic
    [inversemat, stat,reconstructed] = inversemodel(leadfield,'prctile',prc(i));
    source_data=inversemat*preprocessed_eeg;
    EEG_recon=leadfield*source_data;
    for ch=1:32
        R=corrcoef(EEG_recon(ch,:)',preprocessed_eeg(ch,:)');
        corr_all32(i,ch)=R(1,2);
    end
    corr_all(i)=mean(corr_all32(i,:));
    leakage=inversemat*leadfield;
    norm_leakage=leakage./(diag(leakage)*ones(1,size(leakage,2)));
    leak_all(i)=sum(abs(triu(norm_leakage)),'all')-size(norm_leakage,1);
    toc
end

figure
clf
hold on;
yyaxis left; 
plot(prc,corr_all); ylabel({'average corrcoef','Correlation between original and reconstructed EEG'});
xlabel('prctile');
yyaxis right;
plot(prc,log(leak_all)); ylabel({'quantified leakage','log(sum of off-diagonal leakages)'});
grid on;
grid minor;
title('using inversemat.m')

% figure;
% plotx(1:32,corr_all32')
% legend({['prc=1'],['prc=5'],['prc=10'],['prc=25'],['prc=50'],['prc=75']});
% title({['Correlation between original and reconstructed EEG: ']});
% subtitle('using matlab function : inversemodel(leadfield,''prctile'',prc)')
% xlabel('chan'); ylabel('corrcoef');
% grid on;

%% pick prc=1 and pca within regions
load('preprocessed_eeg.mat')
load('20220804_scale_0.05_depth_0.8.mat')
cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat, stat,reconstructed] = inversemodel(leadfield,'prctile',1);
source_data=inversemat*preprocessed_eeg;

tic
[COEFF, SCORE, LATENT] = pca(source_data','Centered',false);
toc

subject_ID='20220804';
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
% bar(source_labels)
% min(unique(source_labels))
fra_eigenvalues=zeros(1,max(unique(source_labels)));
% ave_source_data=zeros(size(source_data,2),max(unique(source_labels)));
% ave_source_coor=zeros(max(unique(source_labels)),3)
ave_source_data=[];
agr_source_data=[];
ave_source_coor=[];ave_source_label=[];
tic
parfor sr=1:max(unique(source_labels))
    I=find(source_labels==sr);
%     figure;plotx(source_data(I,:)')
    if ~isempty(I)
        [COEFF, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
    %     figure;plotx(SCORE(:,1)*COEFF(:,1)');
        fra_eigenvalues(sr)=LATENT(1)/sum(LATENT);
        if fra_eigenvalues(sr) > 0.5
            ave_source_data=[ave_source_data mean(source_data(I,:)',2)];
            agr_source_data=[agr_source_data SCORE(:,1)];
            ave_source_coor=[ave_source_coor; mean(source_fsaverage(I,:),1)];
            ave_source_label=[ave_source_label; sr];
        end
    end
end
toc
figure;
plot(1:max(unique(source_labels)),fra_eigenvalues,'.','MarkerSize',12);
yline(0.5)
xlabel('ROI labels');ylabel('fraction of 1st eigenvalue')
figure;
plot3(ave_source_coor(1:230,1),ave_source_coor(1:230,2),ave_source_coor(1:230,3),'r.');hold on;
plot3(ave_source_coor(231:end,1),ave_source_coor(231:end,2),ave_source_coor(231:end,3),'b.');
%%
sr=68
I=find(source_labels==sr);
  figure;plotx(source_data(I,:)')
if ~isempty(I)
    [COEFF, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
%     figure;plotx(SCORE(:,1)*COEFF(:,1)');
    fra_eigenvalues(sr)=LATENT(1)/sum(LATENT);
    if fra_eigenvalues(sr) > 0.5
        ave_source_data=[ave_source_data mean(source_data(I,:)',2)];
        ave_source_coor=[ave_source_coor; mean(source_fsaverage(I,:),1)];
        ave_source_label=[ave_source_label; sr];
    end
end
%%
% save('agr_source.mat','ave_source_data','agr_source_data','ave_source_coor','ave_source_label');
save('agr_source.mat','agr_source_data','ave_source_coor','ave_source_label');