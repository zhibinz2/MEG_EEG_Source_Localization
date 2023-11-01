function [X_c] = r2c(X_r)
% This function convert the real value covariance to complex value
% using equation 2.21 in the book of Peter J. Schreier, Louis L. Scharf
% -Statistical Signal Processing of Complex-Valued Data_ The Theory of Improper and Noncircular Signals (2010)
n = length(X_r)/2;
Ruu=X_r(1:n,1:n);
Rvv=X_r(n+1:2*n,n+1:2*n);
Ruv=X_r(1:n,n+1:2*n);
Rvu=X_r(n+1:2*n,1:n);
X_c=Ruu+Rvv+1i*(Ruv'-Rvu');
end