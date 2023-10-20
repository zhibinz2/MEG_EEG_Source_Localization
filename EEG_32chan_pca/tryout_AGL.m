%%
clear
cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/
load('./util/dataNeededForSim_MEG.mat') 
powerForWeightedL2 = .5;  % weighting the L2 norm inverse solution
samps = 2400;
SNR = 25;
penUsed = 1; % penalization applied to get inverse of lead field
addpath(genpath('./util'))
addpath(genpath('../AGL/'))
run genCov_CMVN_SC.m
run allInverses_MEG_dip.m 

GforFit =[(eye(length(origNetwork)) +double(origNetwork)), ... 
    double(origNetwork) ; double(origNetwork), (eye(length(origNetwork)) +double(origNetwork))]; % add the dignonal
figure;clf;imagesc(origNetwork);colorbar;colormap('jet')
figure;imagesc(GforFit);colorbar;colormap('jet')


samps = 2400;

allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]); 
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);

% allLambdas = 0.01
% allLambdasOut = 0.6

% simulate data
data =mvnrnd(zeros(length(Q),1),(Q)\eye(size(Q)),samps)';

% store all data:
allData= data;
allOrigPrec= Q;
% generate noise to add to data
tmp = (Q)\eye(size(Q));
noiseAmt = trace(tmp) /(length(tmp)* 10 ^ (SNR/10));
noiseAmtCov = noiseAmt * eye(size(tmp));
noises = mvnrnd(zeros(length(noiseAmtCov),1),(noiseAmtCov),samps);
noisesComp = transpose(noises(:,1:114) +1j*noises(:,115:228));
% create complex valued data to add noise and forward model it to the scalp.
compData = data(1:length(useAreas),:) +1j*data(length(useAreas)+1:2*length(useAreas),:);
% forward model data to the MEG:
chanCompData = L*(compData+noisesComp); % adding complex valued noise and forward modeling
tmp = cat(1, real(chanCompData), imag(chanCompData)); % getting real-valued isomorphism back
tmp = cov(tmp'); % true covariance of the real-valued observation
% # but tmp is not used aterwards

% define source data:
sourceData = squeeze(allOpInv(penUsed,:,:))'*(chanCompData);%+noisesComp

sourceDataReal = cat(1,real(sourceData),imag(sourceData));% 2400x228
sourceDataReal = sourceDataReal*(1/mean(abs(sourceDataReal(:)))); % normalize data
datareshaped = permute(reshape(sourceDataReal', 4,samps/4, size(sourceDataReal,1)),[1,3,2]); % split into 4 ensambles: 4x228x600
% examine the reshape_data
ch=1; n_split=4;sam_len=samps; sam_size=sam_len/n_split; sam_range=1:sam_size; ensam=1;
figure;subplot(121);plot(1:sam_size,sourceDataReal(ch,sam_range),'r');
subplot(122); plot(1:sam_size,squeeze(datareshaped(ensam,ch,sam_range)),'b');

sourceDataReal = cat(1,real(sourceData),imag(sourceData));    
sourceDataReal = [sourceDataReal*(1/mean(abs(sourceDataReal(:))))]'; % normalize data
ch=1; n_split=2; sam_len=samps; sam_size=sam_len/n_split; sam_range=1:sam_size;ensam=1; n_sr=size(sourceDataReal,2);
datareshaped =reshape(sourceDataReal, samps/n_split, n_split, n_sr);
datapermuted = permute(datareshaped,[2,3,1]); % split into 2 ensambles: 2x228x1200
figure;subplot(121);plot(1:sam_size,sourceDataReal(sam_range,ch),'r');
subplot(122); plot(1:sam_size,squeeze(datapermuted(ensam,ch,sam_range)),'b');

tic
% apply adaptive graphical lasso
clear networkPrecCompTrue penInCompTrue penOutCompTrue allDevsReturnTrue
data=datapermuted;
SC=origNetwork;
flagForReal=0;
% [networkPrecCompTrue, penInCompTrue, penOutCompTrue,~,allDevsReturnTrue] = estBestPenalizationQUIC(... 
%         data,SC,allLambdas,allLambdasOut, flagForReal);
% figure;imagesc(networkPrecCompTrue);colorbar;colormap('jet')
% toc 1 m or 8 s
[network,penalizationIn, penalizationOut,minInd,allDevsReturn] = New_estBestPenalizationQUI( ...
    data,SC,allLambdas,allLambdasOut)
toc % # 1 m or 4.2s
figure;imagesc(network);colorbar;colormap('jet')

%% Decifer AGL function -- estBestPenalizationQUIC
[network,penalizationIn,penalizationOut,minInd,allDevsReturn] = New_estBestPenalizationQUI(data,SC,allLambdas,allLambdasOut)
% output: network,penalizationIn, penalizationOut,minInd,allDevsReturn

% input: data,SC,allLambdas,allLambdasOut, flagForReal

% flagForReal=0;
% data=datareshaped;
data=datapermuted;
% SC=origNetwork;
allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]); 
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
min_LamdaIn = min(allLambdas);
n_Lambdas=length(allLambdas);
% origNetwork = double(~eye(length(SC))) .* SC > 0;
% figure;imagesc(ans);colorbar;colormap('jet')
% figure;imagesc(SC);colorbar;colormap('jet')
% figure;imagesc(origNetwork);colorbar;colormap('jet')
origNetwork=SC;

%     if ~flagForReal
        % add the diagonal
        % GforFit =[(eye(length(origNetwork)) +double(origNetwork)),double(origNetwork) ; double(origNetwork), (eye(length(origNetwork)) +double(origNetwork))];
        % without the diagonal
        GforFit =[double(origNetwork),double(origNetwork) ; double(origNetwork), double(origNetwork)];
        % figure;imagesc(GforFit);colorbar;colormap('jet')
%     else
%         GforFit = (eye(length(origNetwork)) +double(origNetwork));
%     end
n_ensam=size(data,1);

allDevs = zeros(length(allLambdas),length(allLambdasOut),n_ensam,n_ensam); % 13x13x4x4 deviance or 13x13x2x2
tic
for mins = 1:n_ensam   % 4 or 2 ensembles for cross validation
    
    parfor lambda = 1:n_Lambdas % loop through lambdas In. Let's pick: lambda = 2
        dataCov = cov(squeeze(data(mins,:,:))'); % 228 x 228 or 888 x 888
        % scalingVal = max(max(triu(abs(dataCov),1)));
        scalingVal = max(triu(abs(dataCov),1),[],'all'); % maximum value of the cov
        tmpDev = zeros(length(allLambdasOut),n_ensam); % 13x4 deviance 
        % cnt = lambda; % indicate which lambda In
        
        for lambdaOut = 1:length(allLambdasOut) % loop through lambdas Out. Let's pick: lambdaOut = 3
            current_lamdaOut=allLambdasOut(lambdaOut);
            current_lamdaIn=allLambdas(lambda);
            if current_lamdaOut < current_lamdaIn % if Out < In, bypass this loop. We want Out > In
                continue
            end
            % offSCedgesLambda = allLambdasOut(lambdaOut); % the current Out lambda value
            % partial coh
            In_network = current_lamdaIn*scalingVal*GforFit;
            Out_network = current_lamdaOut*scalingVal*(double(~(GforFit))-eye(length(GforFit)));
            diag_network = min_LamdaIn*scalingVal*eye(length(GforFit));
            [X] = QUIC('default', dataCov, ... 
                In_network + Out_network + diag_network, ... 
                 1e-4, 0, 200); 


%                 Here we set tol (convergence threshold) = 1e-4
%                 msg (verbosity of print statistics) = 0
%                 maxIter (maximum number of Newton iterations to
%                 execute)=200

%             In_network = allLambdas(lambda)*scalingVal*GforFit;
%             figure;imagesc(In_network);colorbar;colormap('jet');
%             Out_network = offSCedgesLambda*scalingVal*(double(~(GforFit))-eye(length(GforFit)));
%             figure;imagesc(Out_network);colorbar;colormap('jet');
%             figure;imagesc(In_network  +  Out_network);colorbar;colormap('jet');
%             diag_network = min(allLambdas)*scalingVal*eye(length(GforFit));
%             figure;imagesc(diag_network);colorbar;colormap('jet');
%             figure;imagesc(In_network+Out_network+diag_network);colorbar;colormap('jet'); title('L: matrix of the regularization parameters')

            % inverse covariance matrix output
%             figure;imagesc(X);colorbar;colormap('jet');clim([-1 1]);
%             figure;imagesc(abs(X));colorbar;colormap('jet');clim([-1 1]);
%             figure;imagesc(abs(X)>0);colorbar;colormap('jet');clim([0 1]);title('newG');
            newG = abs(X)>0 ;
%             figure;imagesc(newG);colorbar;colormap('jet');clim([0 1]);title('newG');   

%             if ~flagForReal
                newG1 = reduce2nNetwork(newG); % add up the four quadrants
%                 figure;imagesc(newG1);colorbar;colormap('jet');clim([0 1]);title('newG1');   
                % GforFit_new = [newG1, newG1.*double(~eye(length(newG1))); newG1.*double(~eye(length(newG1))), newG1]; % recombine four quadrants
%                 figure;imagesc(GforFit_new);colorbar;colormap('jet');clim([0 1]);title('GforFit_new');
%             else
%                 GforFit_new = newG;
%             end

%             Pweighted = ggmFitHtf(dataCov+ eye(length(dataCov))*min(allLambdas)*max(max(triu(abs(dataCov),1))),GforFit_new);
            % Pweighted = ggmFitHtf(dataCov+ diag_network,GforFit_new); % did not used it
%             figure;imagesc(Pweighted);colorbar;colormap('jet');clim([0 1]);title('Pweighted');

            
            % cross validation
            cnt_min = 1;
            dataCovs = zeros(length(setdiff(1:n_ensam,mins)),size(dataCov,2),size(dataCov,2)); % 3 x 228 x 228
            for mins1 = setdiff(1:n_ensam,mins) % cross validation, loop through 3/4 of the data
                dataCovs(cnt_min,:,:) = cov(squeeze(data(mins1,:,:))'); % 1/4 of the data, righ side 228 x 228
                cnt_min= cnt_min+1;
            end
            
            for  mins1 = setdiff(1:n_ensam,mins)
                % Compute deviance
                % tmpDev(lambdaOut,mins1) = deviance(squeeze(mean(dataCovs,1)), X); 
                S = squeeze(mean(dataCovs,1)); % S = covariance of data (3/4)
                [~,s,~] = svd(X); % X = inverse covariance of model
                tmp = diag(s);
                tmp = tmp(diag(s)>eps);
                logDetD = sum(log(tmp));
                tmpDev(lambdaOut,mins1) = trace(S*X)-logDetD;
            end
            % cnt = cnt + 1; % move on to next lambda In
        end
        allDevs(lambda,:,mins,:) = round(tmpDev,2);
    end
end
toc % 1 m or 32s

allDevs(allDevs==0) = NaN;
% allDevs(imag(allDevs)~=0) = NaN; % if there is complex value (there is none)
% tmp = nanmean(squeeze(nanmean(allDevs,4)),3);
tmp = nanmean(nanmean(allDevs,4),3); % average deviance for each combination of lambda
% tmp1 = nansum(squeeze(nansum(~isnan(allDevs),4)),3);
tmp1 = nansum(nansum(~isnan(allDevs),4),3); % Replaced with the number of cv for each combination of lambda
% check to ensure the penalization chosen has succeeded in creating posDef matrices across all CV
tmp = tmp.*(tmp1==n_ensam*(n_ensam-1)); % correct, nothing changed
tmp(tmp==0) = NaN;
numLocs = n_ensam*(n_ensam-1) - 1; % initiate while loop

% following code ensures we pick some value for the penalization. 
while sum(isnan(tmp(:))) == size(tmp,1)^2 && numLocs ~= 0 % false. skipped
    tmp = nanmean(squeeze(nanmean(allDevs,4)),3);
    tmp1 = nansum(squeeze(nansum(~isnan(allDevs),4)),3);
    tmp = tmp.*(tmp1>numLocs); % correct, nothing changed
    tmp(tmp==0) = NaN; %  nothing changed
    numLocs  = numLocs -1; % update while loop    
end
% nothing changed after this while loop

% find the smallest deviance = smallest penalization
[~,ind] = min(tmp(:)); % linear index
[I1, I2] = ind2sub([size(tmp,1), size(tmp,2)],ind); % convert to subscript index in the 13x13 matrix
allDevsReturn = (tmp); % nothing change

[~,minInd] = min(squeeze(nanmean(allDevs(I1,I2,:,:),3))); % output minInd = index of smallest cv number(not needed)
penalizationIn = allLambdas(I1); % identify the lamda_In value with smallest penalization

data0 = permute(data,[3,1,2]); 
data1 = reshape(data0,size(data0,1)*size(data0,2),size(data0,3)); % combine 4 ensambles back together
% The above step is not needed, if we just use the original sourceDataReal as input
figure;subplot(121);plot(1:1200,data0(:,1,1),'r');
subplot(122);plot(1:1200,data1(1:1200,1),'b');


% compute the covariance of the whole data
dataCov = cov(data1); % 228 x 228
offSCedgesLambda = allLambdasOut(I2); % pick the lamda_Out value with smallest penalization
penalizationOut = allLambdasOut(I2); % same as the the above (used as output)

scalingVal = max(triu(abs(dataCov),1),[],'all'); % maximum value of the cov
In_network = penalizationIn*scalingVal*GforFit;
Out_network = offSCedgesLambda*scalingVal*(double(~(GforFit))-eye(length(GforFit)));
diag_network = min_LamdaIn*scalingVal*eye(length(GforFit));
X = QUIC('default', dataCov, ... 
    In_network +  Out_network + diag_network, ...
     1e-4, 0, 200); % output


newG = abs(X)>0 ; % convert X to bivariate (logical)
% if ~flagForReal
    newG1 = reduce2nNetwork(newG); % add up the four quadrants
%     GforFit_new = [newG1, newG1.*double(~eye(length(newG1))); newG1.*double(~eye(length(newG1))), newG1];
% else
%     GforFit_new = newG;
% end

% network = ggmFitHtf(dataCov+ eye(length(dataCov)) *min(allLambdas)  * max(max(triu(abs(dataCov),1))),GforFit_new);
network = X; 

%% My data
% load source data
clear 
cd /home/zhibinz2/Documents/GitHub/Cleaned_data/cortical_source_data/20220713
load('subj1_tr_1.mat')
cd ../
load('corti_ave_source_labl.mat')
load('corti_ave_source_coor.mat')
source_labels=corti_ave_source_labl{1,1,1};
source_coor=corti_ave_source_coor{1,1,1};

sampl_rate=2000;
srnew = 200;
downsample = 10;
passbands = [1 3; 3.5 6.5; 7 10; 10.5 13.5; 14 20; 21 29; 30 49.5];
bandlabels = {'Delta', 'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2', 'Gamma'};
% Select frequency
freqBand=6;
attenuation=60;
passFreq1 = passbands(freqBand,1);
passFreq2 = passbands(freqBand,2);
d = designfilt('bandpassiir','FilterOrder',20, ...
    'PassbandFrequency1',passFreq1,'PassbandFrequency2',passFreq2, ...
    'StopbandAttenuation1',attenuation,'PassbandRipple',0.2, ...
    'StopbandAttenuation2',attenuation,'SampleRate',srnew);
% test on all channel
% dataIn = double(agr_source_data(151001:153000,:));
dataIn = double(agr_source_data);
downsample_data=resample(double(dataIn),1,downsample,'Dimension',1);
filterd_data = filter(d,downsample_data);
hilbertdata = hilbert(filterd_data'); 
ampdata = abs(hilbertdata); 
angledata = angle(hilbertdata);
% coherence / cross spectra
coh=cov(hilbertdata');
covMat = real2Complex(coh,1);% it doesn't convert to real 
% correlation
cor=corrcoef(ampdata');
% imagesc(cor);colorbar;colormap("jet");title('amplitude correlation');subtitle('beta2')

sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));    
sourceDataReal = sourceDataReal*(1/mean(abs(sourceDataReal(:)))); % normalize data

n_split=2; % number of ensambles
sam_len=floor(size(hilbertdata,2)/n_split); % sample length in each ensamble
sam_size=sam_len*n_split;
n_sr=size(sourceDataReal,1); % number of sources/variables

% remove a few samples at the end
sourceDataReal = [sourceDataReal(:,1:sam_size)]';
datareshaped = reshape(sourceDataReal, sam_len, n_split, n_sr);
datapermuted = permute(datareshaped,[2,3,1]);
% imagesc(corrcoef(sourceDataReal));colorbar;colormap("jet");title('correlation');subtitle('beta2')

% examine the permuted data
ch=1; ensam=1; sam_range=1:sam_len;
figure;subplot(121);plot(sam_range,sourceDataReal(sam_range,ch),'r');
subplot(122); plot(sam_range,squeeze(datapermuted(ensam,ch,sam_range)),'b');

% load SC
load('../../Virtual-Tractography/ForZhibin/processed_data/scale250_Connectome.mat');
SC_tr=logical(fc(source_labels,source_labels));
addpath(genpath('/../../AdaptiveGraphicalLassoforParCoh/Simulations/util'))
addpath(genpath('../../AdaptiveGraphicalLassoforParCoh/AGL'))

allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]); 
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);

% allLambdas = 0.01;
% allLambdasOut = 0.6;

addpath(genpath('/home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util'))
open reduce2nNetwork

data=datapermuted;
SC=SC_tr;

tic
% apply adaptive graphical lasso
clear network networkPrecCompTrue penInCompTrue penOutCompTrue minInd allDevsReturnTrue
% [networkPrecCompTrue, penInCompTrue, penOutCompTrue,~,allDevsReturnTrue] ...
%         = estBestPenalizationQUIC(datareshaped, SC_tr, allLambdas, allLambdasOut, 0);
% toc
% # 11448.750612 seconds.
% figure;imagesc(networkPrecCompTrue);colorbar;colormap('jet')
[network,penalizationIn, penalizationOut,minInd,allDevsReturn] = ...
    New_estBestPenalizationQUI(datapermuted,SC_tr,allLambdas,allLambdasOut);
toc % 1357s -> 1034s 
figure;imagesc(network);colorbar;colormap('jet');
title('AGL output')

%% https://www.mathworks.com/help/parallel-computing/run-matlab-functions-on-a-gpu.html

%% https://www.mathworks.com/help/parallel-computing/run-a-batch-job.html
job = batch('myjob','Pool',20);



