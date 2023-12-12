cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/organize_62stroke
clear
matfile_names=[0:1:60]';
days_poststroke=cat(1,readcell('lesion_site_62subj_sort.xls','Sheet','Sheet1','Range','F2:F12'), ...
    readcell('lesion_site_62subj_sort.xls','Sheet','Sheet1','Range','F14:F63'));
cort_subc=cat(1,readcell('lesion_site_62subj_sort.xls','Sheet','Sheet1','Range','H2:H12'), ...
    readcell('lesion_site_62subj_sort.xls','Sheet','Sheet1','Range','H14:H63'));
FM_scores=cat(1,readcell('lesion_site_62subj_sort.xls','Sheet','Sheet1','Range','J2:J12'), ...
    readcell('lesion_site_62subj_sort.xls','Sheet','Sheet1','Range','J14:J63'));
ages=cat(1,readcell('lesion_site_62subj_sort.xls','Sheet','Sheet1','Range','I2:I12'), ...
    readcell('lesion_site_62subj_sort.xls','Sheet','Sheet1','Range','I14:I63'));
lesion_sizes=cat(1,readcell('lesion_site_62subj_sort.xls','Sheet','Sheet1','Range','K2:K12'), ...
    readcell('lesion_site_62subj_sort.xls','Sheet','Sheet1','Range','K14:K63'));

cor_ind=strcmp('cort',cort_subc);
subc_ind=strcmp('subc',cort_subc);

days_poststroke=cell2mat(days_poststroke);
ages=cell2mat(ages);
lesion_sizes=cell2mat(lesion_sizes);
FM_scores=cell2mat(FM_scores);

save('demographics.mat','matfile_names','days_poststroke','cor_ind',"subc_ind",'FM_scores','lesion_sizes','ages');

%% degree_histogram in python