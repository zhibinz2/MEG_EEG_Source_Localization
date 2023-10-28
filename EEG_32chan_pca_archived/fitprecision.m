function [X] = fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov)
% This function fit precison for each trial
GforFit =[double(SC),double(SC) ; double(SC), double(SC)]; % boolean

% maximum value of the cov
scalingVal = max(triu(abs(dataCov),1),[],'all');

In_network = penalizationIn*scalingVal*GforFit;
Out_network = penalizationOut*scalingVal*(double(~(GforFit))-eye(length(GforFit)));
diag_network = min_LamdaIn*scalingVal*eye(length(GforFit));
X = QUIC('default', dataCov, In_network +  Out_network + diag_network, ...
     1e-4, 0, 200); 

% newG = abs(X)>0 ; % convert X to bivariate (logical)
% network = reduce2nNetwork(newG); % add up the four logical quadrants
end