%% examine original partial correlation 
load Pcoh_all.mat
ses=12;
subj=1;
tr=2;
freq=2;
Pcoh=squeeze(Pcoh_all(ses,subj,tr,freq,:,:));
imagesc(Pcoh);colorbar;colormap('jet');
clim([-0.0001 0.0001]);