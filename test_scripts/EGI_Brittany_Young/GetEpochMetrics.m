
function [OHA,THV,CHV] = GetEpochMetrics(epoch,settings)
X = epoch;
t = size(X,2);
c = size(X,1);
X = X - repmat(nanmean(X,1),c,1);

% overall ratio of timepoints of high amplitude
OHA = nansum(abs(X(:)) > settings.overallThresh)./(t.*c);
% ratio of timepoints of high variance
THV = nansum(bsxfun(@gt, std(X,[],1)', settings.timeThresh), 1) ./t;
% get the number of channels above threshold...
CHV = sum(nanstd(X,[],2) > settings.chanThresh, 1)./c;

end