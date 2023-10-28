function dev = deviance(S,D)
    %S - covariance matrix of data
    %D - inverse covariance matrix of model
%     logDetD = log(det(S*D));
%     if ~isinf(abs(logDetD))
%         dev = - log(det(D));
%     else
        [~,s,~] = svd(D);
        tmp = (diag(s));
        tmp = tmp(tmp>eps);
        logDetD = sum(log(tmp));
        dev = trace(S*D) - logDetD;
%     end
end