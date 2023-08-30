clear
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_streamline_JOHR/
subject_ID='JOHR'; 
% scale=0.01; depth=2;
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
colors9=[red;pink;black;blue;darkgreen;deepyellow;matlab_blue;matlab_orange;matlab_purple];
colors4=[matlab_blue;matlab_orange;deepyellow;darkgreen];

%%
scales={'0.01','0.05','0.1','0.3','0.5'}; depths={'0.8','1','2','4'};
figure;
for s = 1:length(scales)
    subplot(1,length(scales),s)
    legends=cell(1,length(depths));
    hold on;
    for d = 1:length(depths)
        filename=[subject_ID '_scale_' scales{s} '_depth_' depths{d} '.mat'];
        load(filename);
        hold on;
        plot(1:length(corrcoef_diag),corrcoef_diag,'.', 'markersize', 5, 'color',colors4(d,:));
        legends{d}=['depth ' depths{d}];
        clear corrcoef_diag
    end
    hold off;
    legend(legends,'location','southeast')
    title(['scale ' scales{s}])
    xlabel('channel');ylabel('corrcoef');
    clear legends
end
sgtitle([subject_ID ': correlation between original and reconstructed EEG']);
%%
scales={'0.01','0.05','0.1','0.3','0.5'}; depths={'0.8','1','2','4'};
scales3=[0.01 0.05 0.1 0.3 0.5]; depths3=[0.8 1 2 4];
scale_var=nan(1,length(scales)*length(depths));
depth_var=nan(1,length(scales)*length(depths));
corrcoef_var=nan(1,length(scales)*length(depths));
c=1;
for s = 1:length(scales)
    legends=cell(1,length(depths));
    hold on;
    for d = 1:length(depths)
        filename=[subject_ID '_scale_' scales{s} '_depth_' depths{d} '.mat'];
        load(filename);
        corrcoef_var(c)=mean(corrcoef_diag);
        scale_var(c)=scales3(s);depth_var(c)=depths3(d);
        clear corrcoef_diag
        c=c+1;
    end
end

figure;
plot3(scale_var,depth_var,corrcoef_var,'.');
xlabel('scale');ylabel('depth');zlabel('corrcoef-ave');
title(subject_ID);
grid on

