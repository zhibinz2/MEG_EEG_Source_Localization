clear
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_spacing_32chan
subject_ID='20220713'; 
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
%%
scales={'0.01', '0.05', '0.1', '0.3', '0.5'};
depths={'0.8'};
i=1;d=1;
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
%%
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

%% pca
clear
load hald
[COEFF, SCORE, LATENT] = pca(ingredients,'Centered',false);
recon_ingredients= SCORE(:,1:3)*COEFF(:,1:3)';
%% svd
A = [1 2; 3 4; 5 6; 7 8]
[U,S,V] = svd(A)
U*S(:,1)*V(:,1)'
