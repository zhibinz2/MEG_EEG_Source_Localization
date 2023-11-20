clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% cd ../../Cleaned_data/hilbert_chancov
% for ses=1:12
%     mkdir(num2str(seeds(ses,:)));
% end

%% examine all trial
load('chan_prec_all.mat')
ses=12;
subj=2;
tr=12;
freq=5;
imagesc(squeeze(chan_prec_all(ses,subj,tr,freq,:,:)));

for ses=1:12
    for subj=1:2
        for tr=1:12
            for freq=1:5
                imagesc(squeeze(chan_prec_all(ses,subj,tr,freq,:,:)));
                colorbar;
                colormap('jet')
                clim([-1 1]);
                title(['ses ' num2str(ses) ...
                    ' subj ' num2str(subj) ...
                    ' tr ' num2str(tr) ...
                    ' freq ' num2str(freq)])
                pause(0.1);
            end 
        end
    end
end

%% convert back to complex and then normalize to cross spetra
