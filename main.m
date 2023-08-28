%% Loading filtered sono signal
clear, clc
addpath('./RegularizationFunctions/',"./ElastographyFunctions")

baseDir = 'C:\Users\sebas\Documents\MATLAB\Elastography';
%baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\Datasets';
rawPath = [baseDir,'\heterogeneo_sono'];

sonoPath = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\TV\sono';

% Colormap
load('MyColormaps.mat')

%% Subsampling and cropping
% Filtering
swsRange = [2,8];  
move = 'left';

% Cropping
z0 = 4e-3;
zf = 28e-3;

% Looping
Nim = 9; % Number of channels
VibFreqArray = 200:20:360;

fig1 = figure('Position',[100 100 1000 500]);
t1 = tiledlayout(fig1,3,3);
fig2 = figure('Position',[100 100 1000 500]);
t2 = tiledlayout(fig2,3,3);
for iIm = 1:Nim
    % Pre-processing
    load([rawPath, '\Image',num2str(iIm),'\sono.mat']);
    Properties.pitch = 3.0800e-04;
    [sonoFilt,~,~] = process_sono_data(sono,Properties,move,swsRange);

    % Subsampling
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
    
    % Selecting ROI
    zNew = Properties.Depth_S(1:2:end);
    izROI = zNew>z0 & zNew<zf;
    sono = sonoNew(izROI,:,:);
    Properties.Depth_S = zNew(izROI);
    izBmode = Properties.Depth_B>z0 & Properties.Depth_B<zf;
    Properties.Bmode = Properties.Bmode(izBmode,:);
    Properties.Depth_B = Properties.Depth_B(izBmode);
    
    x = Properties.Width_S*1000;
    z = Properties.Depth_S*1000;
    nexttile(t1,iIm);
    imagesc(x,z,sono(:,:,1),1.5*[-1 1])
    colormap(sonomap)
    colorbar
    axis equal
    axis tight
    title(['Sono f_v=',num2str(VibFreqArray(iIm))])

    x = Properties.Width_B*1000;
    z = Properties.Depth_B*1000;
    nexttile(t2,iIm);
    imagesc(x,z,Properties.Bmode)
    colormap gray
    colorbar
    axis equal
    axis tight
    title(['B-mode f_v=',num2str(VibFreqArray(iIm))])
    
    save([sonoPath,'\',num2str(iIm),'.mat'],"sono","Properties");
end

%% --------------------------------------------------------------------- %%
%% 2-DIMENSIONAL TOTAL VARIATION
Nim = 9; % Number of channels
swsPath = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\TV\sws';
ParamsTV.mu = 5;
ParamsTV.tol = 1e-1;
ParamsTV.isotropic = true;

RW_Param.dual = boolean(1); RW_Param.k = 1;
RW_Param.N_window = 20; RW_Param.beta = 1/100000;
RW_Param.tolerance = 1/1000;RW_Param.operator = 'G';
RW_Param.alpha=2;

fig1 = figure('Position',[100 100 1000 500]);
t1 = tiledlayout(fig1,3,3);
fig2 = figure('Position',[100 100 1000 500]);
t2 = tiledlayout(fig2,3,3);
fig3 = figure('Position',[100 100 1000 500]);
t3 = tiledlayout(fig3,3,3);
for iIm = 1:Nim
    load([sonoPath,'\',num2str(iIm),'.mat'])
    fprintf("\nVibration Frequency = %d Hz\n",Properties.VibFreq);
    [swsTV,C] = calculateSWSTV(sono,Properties,ParamsTV);
    fprintf("Number of iterations: %d\n",length(C));
    
    tic
    swsRW = rWave2(sono,Properties,RW_Param);
    t = toc;
    fprintf('Exec. time for R-W: %f\n',t)

    tic
    swsCWT = process_CWT(sono,Properties,[1 16]);
    t = toc;
    fprintf('Exec. time for CWT: %f\n',t)
    
    save([swsPath,'\',num2str(iIm),'.mat'],'swsTV','C','swsRW','swsCWT');

    % Selecting ROI
    x = Properties.Width_S*1000;
    z = Properties.Depth_S*1000;
    SWS_im_range = [2 6];
    
    nexttile(t1,iIm);
    imagesc(x,z,swsTV,SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['SWS from TV, \mu=',num2str(ParamsTV.mu,2),...
        ', f_v=',num2str(Properties.VibFreq)])
    ax = gca; ax.FontSize = 12;

    nexttile(t2,iIm);
    imagesc(x,z,swsRW,SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['SWS from RW, \alpha=',num2str(RW_Param.alpha,2),...
        ', f_v=',num2str(Properties.VibFreq)])
    ax = gca; ax.FontSize = 12;
 
    nexttile(t3,iIm);
    swsCwtIm = squeeze(mean(swsCWT,3));
    imagesc(x,z,swsCwtIm,SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['SWS from CWT, f_v=',num2str(Properties.VibFreq)])
    ax = gca; ax.FontSize = 12;
end

%% Selecting ROI
x0inc = 15; z0 = 11; L = 11; x0back = 1.5;
figure, 
%imagesc(x,z,swsTV,SWS_im_range); colormap turbo
imagesc(Properties.Width_B*1000,Properties.Depth_B*1000,...
    Properties.Bmode); colormap gray
colorbar
axis equal
hold on
rectangle('Position',[x0inc z0 L L], 'LineWidth',2),
rectangle('Position',[x0back z0 L L], 'LineWidth',2, 'EdgeColor','w'),
hold off
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title(['SWS from TV, \mu=',num2str(ParamsTV.mu,2),...
    ', f_v=',num2str(Properties.VibFreq)])
ax = gca; ax.FontSize = 12;

%% Calculating CNR
VibFreqArray = 200:20:360;
[X,Z] = meshgrid(x,z);
maskInc = (X>x0inc & X<x0inc+L & Z>z0 & Z<z0+L);
maskBack = (X>x0back & X<x0back+L & Z>z0 & Z<z0+L);

cnrTV = zeros(1,Nim);
cnrRW = zeros(1,Nim);
cnrCWT = zeros(1,Nim);

meanIncTV = zeros(1,Nim);
meanBackTV = zeros(1,Nim);
stdIncTV = zeros(1,Nim);
stdBackTV = zeros(1,Nim);

meanIncRW = zeros(1,Nim);
meanBackRW = zeros(1,Nim);
stdIncRW = zeros(1,Nim);
stdBackRW = zeros(1,Nim);

meanIncCWT = zeros(1,Nim);
meanBackCWT = zeros(1,Nim);
stdIncCWT = zeros(1,Nim);
stdBackCWT = zeros(1,Nim);
for iIm = 1:Nim
    load([swsPath,'\',num2str(iIm),'.mat']);

    swsInc = swsTV(maskInc);
    swsBack = swsTV(maskBack);
    cnrTV(iIm) = 2*(mean(swsBack) - mean(swsInc))^2 / ...
        (var(swsInc) + var(swsBack));
    meanIncTV(iIm) = mean(swsInc);
    meanBackTV(iIm) = mean(swsBack);
    stdIncTV(iIm) = std(swsInc);
    stdBackTV(iIm) = std(swsBack);  

    swsInc = swsRW(maskInc);
    swsBack = swsRW(maskBack);
    cnrRW(iIm) = 2*(mean(swsBack) - mean(swsInc))^2 / ...
        (var(swsInc) + var(swsBack));
    meanIncRW(iIm) = mean(swsInc);
    meanBackRW(iIm) = mean(swsBack);
    stdIncRW(iIm) = std(swsInc);
    stdBackRW(iIm) = std(swsBack); 

    swsInc = swsCWT(maskInc);
    swsBack = swsCWT(maskBack);
    cnrCWT(iIm) = 2*(mean(swsBack) - mean(swsInc))^2 / ...
        (var(swsInc) + var(swsBack));
    meanIncCWT(iIm) = mean(swsInc);
    meanBackCWT(iIm) = mean(swsBack);
    stdIncCWT(iIm) = std(swsInc);
    stdBackCWT(iIm) = std(swsBack); 
end
%% Comparing CNR
figure('Position', [100 100 400 400]),
plot(VibFreqArray,db(cnrTV),'o-', 'LineWidth',2)
hold on
plot(VibFreqArray,db(cnrRW),'o-', 'LineWidth',2)
plot(VibFreqArray,db(cnrCWT),'o-', 'LineWidth',2)
hold off
legend({'TV','R-WAVE','CWT'}, 'Location','northwest');
ylabel('CNR [dB]'), xlabel('F_v')
grid on
xlim([180 380])

%% Comparing mean SWS and std
figure('Position', [100 100 400 400]),
tiledlayout(1,3)
nexttile
errorbar(VibFreqArray,meanBackTV,stdBackTV, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncTV,stdIncTV, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 380])
title('Total Variation')

nexttile
errorbar(VibFreqArray,meanBackRW,stdBackRW, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncRW,stdIncRW, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 380])
title('R-WAVE')

nexttile
errorbar(VibFreqArray,meanBackCWT,stdBackCWT, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncCWT,stdIncCWT, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 380])
title('CWT')

%% --------------------------------------------------------------------- %%
%% TOTAL NUCLEAR VARIATION
Nim = 9;
for iIm = 1:2:Nim
    % Loading
    load([sonoPath,'\',num2str(iIm),'.mat'])
    [Nz,Nx,~] = size(sono);
    fprintf("\nVibration Frequency = %d Hz\n",Properties.VibFreq);
    
    % Generating matrices
    tic
    [Aim,Bim] = generateSystemTV(sono,Properties);
    t = toc;
    fprintf('Exec. time for generating the system: %f\n',t)
    
    % Re-ordering for col-major
    elemA = reshape(1:Nz*Nx,[Nx,Nz]);
    elemA = elemA';
    icolA = elemA(:);
    Aim = Aim(:,icolA);
    
    % Merging linear systems
    tic
    if iIm == 1
        A = Aim;
        B = Bim;
    else
        [Him,Wim] = size(Aim);
        [Htot,Wtot] = size(A);
        A = [A sparse(Htot,Wim);sparse(Him,Wtot),Aim];
        B = [B;Bim];
    end
    t = toc;
    fprintf('Exec. time for merging the systems: %f\n',t)
end
%% Algorithm
lambda = 0.5;
tau = 0.1;
maxIter = 1000;
tol = 1e-4;
numberEstimators = 5;
stableIter = 50;
[u, cost, error, fide, regul] = pdo_inv_tnv(B, Nz, Nx, A, ...
    lambda, tau, maxIter, tol, numberEstimators, stableIter);
swsTVN = reshape(u,Nz,Nx,numberEstimators);

%% Cost
sws1 = squeeze(swsTVN(:,:,1));
sws2 = squeeze(swsTVN(:,:,5));
x = Properties.Width_S * 1000;
z = Properties.Width_S * 1000;
SWS_im_range = [2,6];

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

%%
subplot(2,2,1)
plot(cost)
title('Cost function')

subplot(2,2,2)
plot(error)
title('Relative error')

subplot(2,2,3)
plot(fide)
title('Fidelity term')

subplot(2,2,4)
plot(regul)
title('Regularization term')


%% --------------------------------------------------------------------- %%
%% USEFUL FUNCTIONS
function [swsTV,C] = calculateSWSTV(sono,Properties,ParamsTV)
% Generates the SWS map from the sonoelastography video
tic
[A,B] = generateSystemTV(sono,Properties);
t = toc;
fprintf('Exec. time for generating the system: %f\n',t)

[Nz,Nx,~] = size(sono);
mask = ones(Nz*Nx,1);
colMajor = false;
tic
[u,C] = IRLS_TV(B,A,ParamsTV.mu,Nz,Nx,ParamsTV.tol,mask,ParamsTV.isotropic,colMajor);
swsTV = reshape(u,Nx,Nz)';
t = toc;
fprintf('Exec. time for solving the system   : %f\n',t)
end

function [A,B] = generateSystemTV(sono,Properties)
% Generates the system of equations of the whole sonoelastography video
[Nz,Nx,~] = size(sono);
A = [];
B = [];
tic
for iz = 1:Nz
    sonoSlice = squeeze(sono(iz,:,:));
    [Az,Bz] = generateWeightMatrix(sonoSlice,Properties);
    [AzDual,BzDual] = generateWeightMatrix(-sonoSlice,Properties);
    Az = [Az;AzDual];
    Bz = [Bz;BzDual];

    if iz == 1
        A = Az;
        B = Bz;
    else
        A = [A sparse(size(A,1),Nx);sparse(size(Az,1),Nx*(iz-1)),Az];
        B = [B;Bz];
    end
end

end