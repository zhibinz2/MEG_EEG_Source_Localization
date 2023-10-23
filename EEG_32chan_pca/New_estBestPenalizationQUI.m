function [X,penalizationIn, penalizationOut,minInd,allDevsReturn] = New_estBestPenalizationQUI(data,SC,allLambdas,allLambdasOut)
% This is a reorganized version of Ani Wodeyar's estBestPenalizationQUIC
% function. It is modified to work with only real value input data.

% this function takes data (ensembles x sources x samples) and estimates a
% model using SC-QUIC for each ensemble and estimates (cross validated) deviance on the
% other ensembles. Using this method it searches for the ideal
% penalizations (lambda1 and lambda2) to apply to the data, and also the ideal epoch to use to
% represent the network represented by the data. It does an upper half
% grid search to identify the ideal penalization (meaning lambda1 >= lambda2).

% For complex data to be used as input in this function,
% real valued matrix data is expected with the real part occupying the first half of sources 
% and imaginary part the second half of sources. 
% See the following lines of code (dataComplex is complex valued, say from fourier transform
% and in the ensembles x sources x samples organization): 
% tmp = real(squeeze(dataComplex(:,:,:)));
% tmp1 = imag(squeeze(dataComplex(:,:,:)));
% tmpCompToRealData = cat(2,tmp,tmp1); 
% For more details refer to Schreir and Scharf, 2010.

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% data: ensembles x sources x samples matrix
% SC: the structure of the penalization, named SC for structural
% connectome. It should be converted to bivariate boolean to be used as
% input.
% allLambdas: the set of lambda values to pass through, will require some trial and error to determine optimal set, assumed to be sorted
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
% network: the precision estimated from the data for the ensemble with
% least deviance in CV and converted to boolean
% penalizationIn: lambda1 penalization
% penalizationOut: lambda2 penalization
% minInd: ensemble selected
% allDevsReturn: the deviance values calculated during cross validation

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zhibin Zhou, Ani Wodeyar
% 10/19/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_LamdaIn = min(allLambdas);
n_ensam=size(data,1); % number of ensambles
n_Lambdas=length(allLambdas); % number of lambda values
GforFit =[double(SC),double(SC) ; double(SC), double(SC)]; % boolean

% initiate the deviance
allDevs = zeros(length(allLambdas),length(allLambdasOut),n_ensam,n_ensam);
for mins = 1:n_ensam   % loop through ensembles for cross validation
    dataCov = cov(squeeze(data(mins,:,:))');
    scalingVal = max(triu(abs(dataCov),1),[],'all'); % maximum value of the cov
    parfor lambda = 1:n_Lambdas % loop through lambdas In.
        tmpDev = zeros(length(allLambdasOut),n_ensam); % 13x4 deviance 
        for lambdaOut = 1:length(allLambdasOut) % loop through lambdas Out.
            current_lamdaOut=allLambdasOut(lambdaOut);
            current_lamdaIn=allLambdas(lambda);
            if current_lamdaOut < current_lamdaIn % if Out < In, bypass this loop. We want Out > In
                continue
            end
            In_network = current_lamdaIn*scalingVal*GforFit;
            Out_network = current_lamdaOut*scalingVal*(double(~(GforFit))-eye(length(GforFit)));
            diag_network = min_LamdaIn*scalingVal*eye(length(GforFit));
            X = QUIC('default', dataCov,In_network + Out_network + diag_network, ... 
                 1e-4, 0, 200);
%             newG = abs(X)>0; % convert to boolean
%             newG1 = reduce2nNetwork(newG); % add up the four logical quadrants

            % cross validation
            cnt_min = 1;
            dataCovs = zeros(length(setdiff(1:n_ensam,mins)),size(dataCov,2),size(dataCov,2)); % 3 x 228 x 228
            for mins1 = setdiff(1:n_ensam,mins) % cross validation, loop through 3/4 of the data
                dataCovs(cnt_min,:,:) = cov(squeeze(data(mins1,:,:))'); % 1/4 of the data, righ side 228 x 228
                cnt_min= cnt_min+1;
            end

            for  mins1 = setdiff(1:n_ensam,mins)
                % Compute deviance
                S = squeeze(mean(dataCovs,1)); % S = covariance of data (3/4)
                [~,s,~] = svd(X); % X = inverse covariance of model
                tmp = diag(s);
                tmp = tmp(diag(s)>eps);
                logDetD = sum(log(tmp));
                tmpDev(lambdaOut,mins1) = trace(S*X)-logDetD;
            end 
        end
        allDevs(lambda,:,mins,:) = round(tmpDev,2);
    end
    clear dataCov
end

allDevs(allDevs==0) = NaN;

% average deviance for each combination of lambda
allDevsReturn = nanmean(nanmean(allDevs,4),3); % average deviance for each combination of lambda

% find the smallest deviance = smallest penalization
[~,~] = min(allDevsReturn(:)); % linear index
[I1, I2] = ind2sub([size(allDevsReturn,1), size(allDevsReturn,2)],ind); % convert to subscript index in the 13x13 matrix

% identify the lamda_In value with minimum penalization
penalizationIn = allLambdas(I1); 
% pick the lamda_Out value with minimum penalization
penalizationOut = allLambdasOut(I2); 
% output minInd = index of which emsamble having the minimum penalization
[~,minInd] = min(squeeze(nanmean(allDevs(I1,I2,:,:),3))); 

% combine the ensambles back together
data0 = permute(data,[3,1,2]); 
data1 = reshape(data0,size(data0,1)*size(data0,2),size(data0,3)); 
clear data0

% compute the covariance of the whole data
dataCov = cov(data1); % 228 x 228
% maximum value of the cov
scalingVal = max(triu(abs(dataCov),1),[],'all');

In_network = penalizationIn*scalingVal*GforFit;
Out_network = penalizationOut*scalingVal*(double(~(GforFit))-eye(length(GforFit)));
diag_network = min_LamdaIn*scalingVal*eye(length(GforFit));
X = QUIC('default', dataCov, In_network +  Out_network + diag_network, ...
     1e-4, 0, 200); 
% newG = abs(X)>0 ; % convert X to bivariate (logical)
% network = reduce2nNetwork(newG); % add up the four logical quadrants
