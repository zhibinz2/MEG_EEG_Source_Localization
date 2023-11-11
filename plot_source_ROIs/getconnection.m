function [connected_rois] = getconnection(fc,ROIs)
% This function gets the labels of the ROIs connected to the input ROIs
% Inputs: 
% fc - full connectome
% ROIs - the interested ROI labels for which we want to know what other ROIs are connected to them
% Outputs:
% connected_rois - the ROIs are connected to the input ROIs

fc_bool=logical(fc);

connected_ones=zeros(size(fc));
connected_ones(ROIs,:)=1;
connected_ones(:,ROIs)=1;
% imagesc(connected_ones)
% length(find(connected_ones));
non_connected_ones=~connected_ones;
% imagesc(non_motor_ones)
% imagesc(fc_bool)
fc_bool(find(non_connected_ones))=0;
% imagesc(fc_bool)

[r,c,v]=find(fc_bool); % r and c are the index of rows and colums, v is the values of ones

% all connected ROIs labels
connected_rois=unique([r;c]);

% plot all sources connected
% [bool_a,ind_b] = ismember(source_labels,connected_rois);

end

