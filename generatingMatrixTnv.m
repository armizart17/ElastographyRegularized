%% Loading filtered sono signal
clear, clc
addpath('./RegularizationFunctions/',"./ElastographyFunctions")
baseDir = 'C:\Users\sebas\Documents\MATLAB\Elastography';
%baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\Datasets';
sonoPath = [baseDir,'\heterogeneo_sono'];
swsRange = [2,8];  
move = 'left';

% Colormap
load('MyColormaps.mat')

% Pre-processing
load([sonoPath, '\Image8\sono.mat']);
Properties.pitch = 3.0800e-04;
[sonoFilt,~,~] = process_sono_data(sono,Properties,move,swsRange);
%% Subsampling in depth
R = 2; % DECIMATION FACTOR
[Nz,Nx,Nt] = size(sonoFilt);
sonoNew = zeros([ceil(Nz/R),Nx,Nt]);
[b,a] = butter(4,1/R);
for ix = 1:Nx
    for it = 1:Nt
        signal = filtfilt(b,a,sonoFilt(:,ix,it));
        sonoNew(:,ix,it) = signal(1:R:end);
    end
end
clear signal
%% Selecting ROI and cleaning workspace
x = Properties.Width_S * 1000;
z = Properties.Depth_S * 1000;

zNew = z(1:2:end);
izROI = zNew>4 & zNew<28;
sonoSub = sonoNew(izROI,:,:);
z = zNew(izROI);
clear sono sonoNew sonoFilt
%% Plotting sono video
figure,
for it = 1:5
    imagesc(x,z,sonoSub(:,:,it),1.5*[-1 1])
    colormap(sonomap)
    colorbar
    axis equal
    axis tight
    pause(1/Properties.FrameRate)
end

%% Generating A matrix
RW_Param.dual = true;
[Nz,Nx,Nt] = size(sonoSub);
A1 = [];
B1 = [];
tic
for iz = 1:Nz
    sonoSlice = squeeze(sonoSub(iz,:,:));
    [Az,Bz] = generateWeightMatrix(sonoSlice,Properties);
    if RW_Param.dual
        [AzDual,BzDual] = generateWeightMatrix(-sonoSlice,Properties);
        Az = [Az;AzDual];
        Bz = [Bz;BzDual];
    end
    if iz == 1
        A1 = Az;
        B1 = Bz;
    else
        A1 = [A1 sparse(size(A1,1),Nx);sparse(size(Az,1),Nx*(iz-1)),Az];
        B1 = [B1;Bz];
    end
end
toc

%% Reordering
tic
elemA = reshape(1:Nz*Nx,[Nx,Nz]);
elemA = elemA';
icolA = elemA(:);
A1 = A1(:,icolA);
toc
clear Az Bz
%% Results with ordered matrix
SWS_im_range = [2,5.5];
mu = 5;
[M,N,~] = size(sonoSub);
tol = 1e-3;
mask = ones(M*N,1);
isotropic = true;
colMajor = true;

[u,G] = IRLS_TV(B1,A1,mu,M,N,tol,mask,isotropic,colMajor);
swsTV = reshape(u,M,N);
%[u,G] = IRLS_TV(B,A,mu,N,M,tol,mask,isotropic,colMajor);
%swsTV = reshape(u,N,M)';

figure('Position', [200 200 400 500]);
subplot(211),
imagesc(x,z,swsTV,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title(['SWS from TV, \mu=',num2str(mu,2)])
ax = gca; ax.FontSize = 12;

subplot(212)
plot(G)
xlabel("# of iterations")
ylabel("Error")
title('Convergence of error')
grid on
axis tight
ax = gca; ax.FontSize = 12;

%% Repeating the process with another image
% Pre-processing
load([sonoPath, '\Image9\sono.mat']);
Properties.pitch = 3.0800e-04;
[sonoFilt,~,~] = process_sono_data(sono,Properties,move,swsRange);

% Subsampling in depth
R = 2; % DECIMATION FACTOR
[Nz,Nx,Nt] = size(sonoFilt);
sonoNew = zeros([ceil(Nz/R),Nx,Nt]);
[b,a] = butter(4,1/R);
for ix = 1:Nx
    for it = 1:Nt
        signal = filtfilt(b,a,sonoFilt(:,ix,it));
        sonoNew(:,ix,it) = signal(1:R:end);
    end
end
clear signal

% Selecting ROI
x = Properties.Width_S * 1000;
z = Properties.Depth_S * 1000;
zNew = z(1:2:end);
izROI = zNew>4 & zNew<28;
sonoSub = sonoNew(izROI,:,:);
z = zNew(izROI);
clear sono sonoNew sonoFilt

%%
figure,
for it = 1:5
    imagesc(x,z,sonoSub(:,:,it),1.5*[-1 1])
    colormap(sonomap)
    colorbar
    axis equal
    axis tight
    pause(1/Properties.FrameRate)
end

%% Generating A matrix
RW_Param.dual = true;
[Nz,Nx,Nt] = size(sonoSub);
A2 = [];
B2 = [];
tic
for iz = 1:Nz
    sonoSlice = squeeze(sonoSub(iz,:,:));
    [Az,Bz] = generateWeightMatrix(sonoSlice,Properties);
    if RW_Param.dual
        [AzDual,BzDual] = generateWeightMatrix(-sonoSlice,Properties);
        Az = [Az;AzDual];
        Bz = [Bz;BzDual];
    end
    if iz == 1
        A2 = Az;
        B2 = Bz;
    else
        A2 = [A2 sparse(size(A2,1),Nx);sparse(size(Az,1),Nx*(iz-1)),Az];
        B2 = [B2;Bz];
    end
end
toc

% Reordering
tic
elemA = reshape(1:Nz*Nx,[Nx,Nz]);
elemA = elemA';
icolA = elemA(:);
A2 = A2(:,icolA);
toc
clear Az Bz
%% Results with ordered matrix
SWS_im_range = [2,5.5];
mu = 5;
[M,N,~] = size(sonoSub);
tol = 1e-3;
mask = ones(M*N,1);
isotropic = true;
colMajor = true;

[u,G] = IRLS_TV(B2,A2,mu,M,N,tol,mask,isotropic,colMajor);
swsTV = reshape(u,M,N);

figure('Position', [200 200 400 500]);
subplot(211),
imagesc(x,z,swsTV,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title(['SWS from TV, \mu=',num2str(mu,2)])
ax = gca; ax.FontSize = 12;

subplot(212)
plot(G)
xlabel("# of iterations")
ylabel("Error")
title('Convergence of error')
grid on
axis tight
ax = gca; ax.FontSize = 12;

%% Joining two channels
[H1,W1] = size(A1);
[H2,W2] = size(A2);
tic
A = [A1 sparse(H1,W2);sparse(H2,W1),A2];
B = [B1;B2];
toc
%% TNV
lambda = 0.5;
tau = 0.1;
maxIter = 1000;
tol = 1e-4;
numberEstimators = 2;
stableIter = 50;
[u, cost, error, fide, regul] = pdo_inv_tnv(B, Nz, Nx, A, lambda, tau, maxIter, tol, numberEstimators, stableIter);

%%
swsTVN = reshape(u,M,N,numberEstimators);
sws1 = squeeze(swsTVN(:,:,1));
sws2 = squeeze(swsTVN(:,:,2));
x = Properties.Width_S * 1000;

figure('Position',[100 100 800 500])
subplot(2,2,1)
imagesc(x,z,sws1,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title(['SWS from TVN, \lambda=',num2str(lambda,2)])
ax = gca; ax.FontSize = 12;

subplot(2,2,3)
imagesc(x,z,sws2,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title(['SWS from TVN, \lambda=',num2str(lambda,2)])
ax = gca; ax.FontSize = 12;

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
