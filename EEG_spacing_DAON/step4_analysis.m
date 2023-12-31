clear
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_spacing_DAON/
subject_ID='DAON'; 
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
colors4=[matlab_blue;matlab_orange;darkgreen];

%%
scales={'0.01','0.05','0.1','0.3','0.5'}; depths={'0.8','1','2','4'};
s=2;d=1;
icos={'2','3','4'}; % icos={'2','3','4','5'}
scale=0.05; depth=0.8;
figure;
legends=cell(1,length(icos));
for i = 1:length(icos)
    hold on;
    filename=[subject_ID '_ico_' icos{i} '_scale_' scales{s} '_depth_' depths{d} '.mat'];
    load(filename);
    hold on;
    plot(1:length(corrcoef_diag),corrcoef_diag,'.', 'markersize', 5, 'color',colors4(i,:));
    legends{i}=['ico ' icos{i}];
end
clear corrcoef_diag
legend(legends,'location','southeast')
xlabel('channel');ylabel('corrcoef');
title([subject_ID ': correlation between original and reconstructed EEG']);
%%
scales={'0.01','0.05','0.1','0.3','0.5'}; depths={'0.8','1','2','4'};
scales3=[0.01 0.05 0.1 0.3 0.5]; depths3=[0.8 1 2 4];
% icos={'2','3','4','5'}; icos3=[2 3 4 5];
icos={'2','3','4'}; icos3=[2 3 4];
% scale_var=nan(1,length(scales)*length(depths));
% depth_var=nan(1,length(scales)*length(depths));
corrcoef_var=nan(1,length(icos3));
c=1;
for i = 1:length(icos)
    legends=cell(1,length(icos));
    hold on;
    filename=[subject_ID '_ico_' icos{i} '_scale_' scales{2} '_depth_' depths{1} '.mat'];
    load(filename);
    corrcoef_var(c)=mean(corrcoef_diag);
%     scale_var(c)=scales3(s);depth_var(c)=depths3(d);
    clear corrcoef_diag
    c=c+1;
end

figure;
scatter(icos3,corrcoef_var);
xlabel('ico');ylabel('corrcoef-ave');
title([subject_ID ': recursively subdivided icosahedron (ico) spacing']);
xlim([1 5]); ylim([0.75 0.78]);
grid on

