%%
clear,clc
addpath("../RegularizationFunctions")
%% Reading image
imBig = imread("papagayo.jpeg");
im = imresize(imBig,0.3);
[H,W,C] = size(im);
imshow(im)
%% Adding noise
y = double(im(:)) + 30*randn(H*W*C,1);
imNoisy = uint8(reshape(y,[H,W,C]));
imshow(imNoisy)
%% Intento de TNV
A = speye(H*W*C);
lambda = 20;
tau = 0.1;
maxIter = 1000;
tol = 1e-4;
numberEstimators = C;
stableIter = 50;
[x, cost, error, fide, regul] = pdo_inv_tnv(y, H, W, A, lambda, tau, maxIter, tol, numberEstimators, stableIter);

%%
imReconstructed = uint8(x);
figure('Position',[100 100 800 500])
subplot(1,2,1)
imshow(imReconstructed)
title('Reconstructed image')

subplot(4,2,2)
plot(cost)
title('Cost function')

subplot(4,2,4)
plot(error)
title('Relative error')

subplot(4,2,6)
plot(fide)
title('Fidelity term')

subplot(4,2,8)
plot(regul)
title('Regularization term')
