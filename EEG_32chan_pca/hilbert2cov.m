%% loop through all cortical source data to compute hilbert and cov in 5 freq
clear

cd /home/zhibinz2/Documents/GitHub/Cleaned_data/cortical_source_data
load('corti_ave_source_labl.mat')
load('corti_ave_source_coor.mat')
source_labels=corti_ave_source_labl{1,1,1};
source_coor=corti_ave_source_coor{1,1,1};

sampl_rate=2000;
srnew = 200;
downsample = 10;
% passbands = [1 3; 3.5 6.5; 7 10; 10.5 13.5; 14 20; 21 29; 30 49.5];
% bandlabels = {'Delta', 'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2', 'Gamma'};
passbands = [3.5 6.5; 7 10; 10.5 13.5; 14 20; 21 29];
bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};
attenuation=60;
filt_ds=cell(1,length(passbands));
for freqBand=1:5
    % Select frequency
    passFreq1 = passbands(freqBand,1);
    passFreq2 = passbands(freqBand,2);
    filt_d = designfilt('bandpassiir','FilterOrder',20, ...
    'PassbandFrequency1',passFreq1,'PassbandFrequency2',passFreq2, ...
    'StopbandAttenuation1',attenuation,'PassbandRipple',0.2, ...
    'StopbandAttenuation2',attenuation,'SampleRate',srnew);
    filt_ds{freqBand}=filt_d;
end

seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

for ses=1:12
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

                % ampdata = abs(hilbertdata); 
                % angledata = angle(hilbertdata);

                sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));
                sourceDataReal = sourceDataReal*(1/mean(abs(sourceDataReal(:)))); % normalize data

                hilbert_dataCov{freq} = cov(sourceDataReal');

                save(['./hilbert_datacov/' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat'],'hilbert_dataCov');
            end
        end
    end
    toc
end
 
                

