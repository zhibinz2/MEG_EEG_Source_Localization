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
