cd /home/zhibinz2/Documents/GitHub/Cleaned_data/hilbert_datacov
load('hilbert_dataCov_all.mat');

%%
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/
for ses=12%12
    tic
    for subj=1:2
        for tr=1:12
            load(['./cortical_source_data/' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat']);

            hilbert_dataCov=cell(1,5);
            for freq=1:5
                
                % resample
                downsample_data=resample(double(agr_source_data),1,downsample,'Dimension',1);
                filterd_data = filter(filt_ds{freq},downsample_data);
                hilbertdata = hilbert(filterd_data'); 

                sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));
                sourceDataReal = sourceDataReal*(1/mean(abs(sourceDataReal(:)))); % normalize data

                hilbert_dataCov{freq} = cov(sourceDataReal');

                save(['./hilbert_datacov/' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat'],'hilbert_dataCov');
            end
        end
    end
    toc
end

%%
clear
tic
cd /ssd/zhibin/Cleaned_sourcedata/cortical_source_data/python_lasso
prec_all_ses_1_12=zeros(12,2,12,5,896,896);

clear prec_all
load('prec_all_ses_0_1.mat')
prec_all_ses_1_2=prec_all;
ses=2;subj=2;tr=5;freq=1;
prec=squeeze(prec_all_ses_1_2(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(1:2,:,:,:,:,:)=prec_all_ses_1_2(1:2,:,:,:,:,:);
clear prec_all_ses_1_2 prec_all

load('prec_all_ses_2_3.mat') % hnlb
prec_all_ses_3_4=prec_all;
ses=4;subj=2;tr=6;freq=5;
prec=squeeze(prec_all_ses_3_4(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(3:4,:,:,:,:,:)=prec_all_ses_3_4(3:4,:,:,:,:,:);
clear prec_all_ses_3_4 prec_all

load('prec_all_ses_4.mat') % ramesh
prec_all_ses_5=prec_all;
ses=10;subj=1;tr=6;freq=2;
prec=squeeze(prec_all_ses_5(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(5,:,:,:,:,:)=prec_all_ses_5(5,:,:,:,:,:);
clear prec_all_ses_5 prec_all

load('prec_all_ses_5.mat') % hnla
prec_all_ses_6=prec_all;
ses=6;subj=2;tr=6;freq=1;
prec=squeeze(prec_all_ses_6(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(6,:,:,:,:,:)=prec_all_ses_6(6,:,:,:,:,:);
clear prec_all_ses_6 prec_all

load('prec_all_ses_6_7.mat');
prec_all_ses_7_8=prec_all;
ses=7;subj=1;tr=5;freq=3;
prec=squeeze(prec_all_ses_7_8(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(7:8,:,:,:,:,:)=prec_all_ses_7_8(7:8,:,:,:,:,:);
clear prec_all_ses_7_8 prec_all

load('prec_all_ses_8.mat') 
prec_all_ses_9=prec_all;
ses=9;subj=1;tr=6;freq=2;
prec=squeeze(prec_all_ses_9(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(9,:,:,:,:,:)=prec_all_ses_9(9,:,:,:,:,:);
clear prec_all_ses_9 prec_all

load('prec_all_ses_9.mat') % hnla
prec_all_ses_10=prec_all;
ses=10;subj=1;tr=6;freq=2;
prec=squeeze(prec_all_ses_10(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(10,:,:,:,:,:)=prec_all_ses_10(10,:,:,:,:,:);
clear prec_all_ses_10 prec_all

load('prec_all_ses_10.mat') % hnla works % ramesh?
prec_all_ses_11=prec_all;
ses=11;subj=2;tr=3;freq=3;
prec=squeeze(prec_all_ses_11(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(11,:,:,:,:,:)=prec_all_ses_11(11,:,:,:,:,:);
clear prec_all_ses_11 prec_all

load('prec_all_ses_11.mat') % ramesh
prec_all_ses_12=prec_all;
ses=12;subj=2;tr=5;freq=2;
prec=squeeze(prec_all_ses_12(ses,subj,tr,freq,:,:));
imagesc(logical(prec));colorbar
prec_all_ses_1_12(12,:,:,:,:,:)=prec_all_ses_12(12,:,:,:,:,:);
clear prec_all_ses_12 prec_all

cd /ssd/zhibin/Cleaned_sourcedata/cortical_source_data/python_lasso
save('prec_all_ses_1_12.mat','prec_all_ses_1_12','-v7.3');
toc 

figure
sumpre_all=nan(12,2,12,5);
for ses=1:12 
    for subj=1:2
        for tr=1:12
            for freq=1:5
                prec=squeeze(prec_all_ses_1_12(ses,subj,tr,freq,:,:));
                imagesc(logical(prec));colorbar;
                sumprec=sum(triu(prec,1),"all");
                title(['sum of prec = ' num2str(sumprec)]);
                sumpre_all(ses,subj,tr,freq)=sumprec;
                subtitle(['ses ' num2str(ses) ' subj ' num2str(subj) ' tr ' num2str(tr) ...
                    ' freq ' num2str(freq)])
                pause(0.1)
            end
        end
    end
end





