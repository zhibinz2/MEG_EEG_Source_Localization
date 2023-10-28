tic
% apply adaptive graphical lasso
clear networkPrecCompTrue penInCompTrue penOutCompTrue minInd allDevsReturnTrue
[network,penalizationIn, penalizationOut,minInd,allDevsReturn] = ...
    New_estBestPenalizationQUI(datapermuted,SC_tr,allLambdas,allLambdasOut)
toc 
figure;imagesc(network);colorbar;colormap('jet')