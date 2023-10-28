function [penalizationIn,penalizationOut,minDev] = penaltyselection(SC,allLambdas,allLambdasOut,dataCovs_op)
% This function runs cross validation for penalty selection

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% dataCov: n x 894 x 894

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
% X: the precision estimated from the data for the ensemble with
% least deviance in CV and converted to boolean
% penalizationIn: lambda1 penalization
% penalizationOut: lambda2 penalization
% minDev: minimum deviance of the Model

%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_LamdaIn = min(allLambdas);
n_Lambdas=length(allLambdas); % number of lambda values

n_ensam=size(dataCovs_op,1); % number of ensambles
GforFit =[double(SC),double(SC) ; double(SC), double(SC)]; % boolean


% initiate the deviance
allDevs = zeros(length(allLambdas),length(allLambdasOut),n_ensam,n_ensam);
for mins = 1:n_ensam   % loop through ensembles for cross validation
    dataCov = squeeze(dataCovs_op(mins,:,:));
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
            X = QUIC('default', dataCov, In_network + Out_network + diag_network, ... 
                 1e-4, 0, 200);
%             newG = abs(X)>0; % convert to boolean
%             newG1 = reduce2nNetwork(newG); % add up the four logical quadrants

            % cross validation
            cnt_min = 1;
            dataCovs = zeros(length(setdiff(1:n_ensam,mins)),size(dataCov,2),size(dataCov,2)); % 3 x 228 x 228
            for mins1 = setdiff(1:n_ensam,mins) % cross validation, loop through 3/4 of the data
                dataCovs(cnt_min,:,:) = squeeze(dataCovs_op(mins1,:,:)); % 1/4 of the data, righ side 228 x 228
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
[minDev,ind] = min(allDevsReturn(:)); % linear index
[I1, I2] = ind2sub([size(allDevsReturn,1), size(allDevsReturn,2)],ind); % convert to subscript index in the 13x13 matrix

% identify the lamda_In value with minimum penalization
penalizationIn = allLambdas(I1); 
% pick the lamda_Out value with minimum penalization
penalizationOut = allLambdasOut(I2); 
% output minInd = index of which emsamble having the minimum penalization
% [~,minInd] = min(squeeze(nanmean(allDevs(I1,I2,:,:),3))); 

end