clear;
tic
A = gpuArray([1 0 1; -1 -2 0; 0 1 -1]);
e = eig(A);
toc

gpuDeviceTable


gpuDevice