%% Loading filtered sono signal
clear, clc
addpath('./RegularizationFunctions/',"./ElastographyFunctions")

baseDir = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\13_09_23';
sonoPath = [baseDir,'\sono'];
croppedPath = [baseDir,'\cropped'];
iqPath = [baseDir,'\iq'];
swsPath = [baseDir,'\sws'];

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
Nim = 6; % Number of channels

fig1 = figure('Position',[100 100 1000 500]);
t1 = tiledlayout(fig1,2,3);
fig2 = figure('Position',[100 100 1000 500]);
t2 = tiledlayout(fig2,2,3);
for iIm = 1:Nim
    % Pre-processing
    load([sonoPath, '\',num2str(iIm),'.mat']);
    load([iqPath, '\',num2str(iIm),'.mat']);
    Bmode = db(IQ(:,:,1));
    clear IQ;
    
    dinf.FRsono = 20;
    Prop.FrameRate = 20; Prop.pitch = dinf.dx;
    Prop.VibFreq = dinf.vibf; Prop.VibFreqOffset = 2;
    move = 'right'; swsRange = [2,7];
    [sonoFilt,~,~] = process_sono_data(sono,Prop,move,swsRange);

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
    zNew = dinf.z(1:2:end);
    izROI = zNew>z0 & zNew<zf;
    valid = sono(:,:,1) ~= 0;
    validNew = valid(1:2:end,:);
    dinf.valid = validNew(izROI,:,:);
    sono = sonoNew(izROI,:,:);
    dinf.z = zNew(izROI);
    
    x = dinf.x*1000;
    z = dinf.z*1000;
    nexttile(t1,iIm);
    imagesc(x,z,sono(:,:,1),1.5*[-1 1])
    colormap(sonomap)
    colorbar
    axis equal
    axis tight
    title(['Sono f_v=',num2str(dinf.vibf)])

    nexttile(t2,iIm);
    imagesc(x,z,Bmode,[110 150])
    colormap gray
    colorbar
    axis equal
    axis tight
    title(['B-mode f_v=',num2str(dinf.vibf)])
    
    save([croppedPath,'\',num2str(iIm),'.mat'],"sono","dinf");
end

%% --------------------------------------------------------------------- %%
%% 2-DIMENSIONAL TOTAL VARIATION
Nim = 6; % Number of channels
ParamsTV.mu = 5;
ParamsTV.tol = 1e-1;
ParamsTV.isotropic = true;

RW_Param.dual = boolean(1); RW_Param.k = 1;
RW_Param.N_window = 20; RW_Param.beta = 1/100000;
RW_Param.tolerance = 1/1000;RW_Param.operator = 'G';
RW_Param.alpha=2;

fig1 = figure('Position',[100 100 1000 500]);
t1 = tiledlayout(fig1,2,3);
fig2 = figure('Position',[100 100 1000 500]);
t2 = tiledlayout(fig2,2,3);
fig3 = figure('Position',[100 100 1000 500]);
t3 = tiledlayout(fig3,2,3);
for iIm = 1:Nim
    load([croppedPath,'\',num2str(iIm),'.mat'])

    fprintf("\nVibration Frequency = %d Hz\n",dinf.vibf);
    Properties.pitch = dinf.dx;
    Properties.VibFreq = dinf.vibf;
    [swsTV,C] = calculateSWSTV(sono,Properties,ParamsTV);
    fprintf("Number of iterations: %d\n",length(C));
    % 
    tic
    swsRW = rWave2(sono,Properties,RW_Param);
    t = toc;
    fprintf('Exec. time for R-W: %f\n',t)

    tic
    swsCWT = process_CWT(sono,Properties,[1 16]);
    t = toc;
    fprintf('Exec. time for CWT: %f\n',t)
    
    save([swsPath,'\',num2str(iIm),'.mat'],...
        'swsTV','C','swsRW','swsCWT','dinf');

    % Selecting ROI
    x = dinf.x*1000;
    z = dinf.z*1000;
    SWS_im_range = [2 6];
    
    nexttile(t1,iIm);
    im = imagesc(x,z,swsTV,SWS_im_range);
    im.AlphaData = dinf.valid;
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['SWS from TV, \mu=',num2str(ParamsTV.mu,2),...
        ', f_v=',num2str(Properties.VibFreq)])
    ax = gca; ax.FontSize = 12;

    nexttile(t2,iIm);
    im = imagesc(x,z,swsRW,SWS_im_range);
    im.AlphaData = dinf.valid;
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
    im = imagesc(x,z,swsCwtIm,SWS_im_range);
    im.AlphaData = dinf.valid;
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['SWS from CWT, f_v=',num2str(Properties.VibFreq)])
    ax = gca; ax.FontSize = 12;
end

%% Selecting ROI
x0inc = -4; z0 = 11; L = 8; x0back = -15;
figure, 
%imagesc(x,z,swsTV,SWS_im_range); colormap turbo
imagesc(dinf.x*1000,dinf.z*1000, Bmode, [110,150]); colormap gray
colorbar
axis equal
hold on
rectangle('Position',[x0inc z0 L L], 'LineWidth',2),
rectangle('Position',[x0back z0 L L], 'LineWidth',2, 'EdgeColor','w'),
hold off
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title('Bmode')
ax = gca; ax.FontSize = 12;

%% Calculating CNR
VibFreqArray = 200:50:450;
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
xlim([180 470])

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
xlim([180 470]), ylim([2 4])
title('Total Variation')

nexttile
errorbar(VibFreqArray,meanBackRW,stdBackRW, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncRW,stdIncRW, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 470]), ylim([2 4])
title('R-WAVE')

nexttile
errorbar(VibFreqArray,meanBackCWT,stdBackCWT, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncCWT,stdIncCWT, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 470]), ylim([2 4])
title('CWT')

%% --------------------------------------------------------------------- %%
%% TOTAL NUCLEAR VARIATION
%selectedImages = [3,4,5,7,8,9];
selectedImages = 1:6;
Nim = length(selectedImages);
for iIm = 1:Nim
    % Loading
    load([croppedPath,'\',num2str(selectedImages(iIm)),'.mat'])
    [Nz,Nx,~] = size(sono);
    Properties.VibFreq = dinf.vibf;
    Properties.pitch = dinf.dx;
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
% larger tau increases num of iterations and makes the results less
% dependendent on lambda
% smaller tau increases agreement between channels
lambda = 1; % +-50%
tau = 0.2;   % +-50%
maxIter = 500;
tol = 3e-4;
stableIter = 50;
Nim = length(selectedImages);
tic
[u, cost, error, fide, regul] = pdo_inv_tnv(B, Nz, Nx, A, ...
    lambda, tau, maxIter, tol, Nim, stableIter);
swsTNV = reshape(u,Nz,Nx,Nim);
toc
save([swsPath,'\TNV2.mat'],'swsTNV','dinf');

%% Results
load([swsPath,'\TNV2.mat'],'swsTNV');
x = dinf.x * 1000;
z = dinf.z * 1000;
SWS_im_range = [2,6];
VibFreqArray = 200:20:360;

figure('Position',[100 100 1200 600]);
%t = tiledlayout(3,3);
t = tiledlayout(2,3);
for iIm = 1:Nim
    nexttile
    im = imagesc(x,z,squeeze(swsTNV(:,:,iIm)),SWS_im_range);
    im.AlphaData = dinf.valid;
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['F_v = ',num2str(VibFreqArray(selectedImages(iIm))),' Hz'])
    ax = gca; ax.FontSize = 12;
end
title(t,['SWS from TNV, \lambda=',num2str(lambda,2),...
    ', \tau=',num2str(tau,2)])

figure,
t = tiledlayout(2,2);
nexttile, plot(cost, 'LineWidth',1.5)
title('Cost function')
axis tight
nexttile, plot(error, 'LineWidth',1.5)
title('Relative error')
axis tight
nexttile, plot(fide, 'LineWidth',1.5)
title('Fidelity term')
axis tight
nexttile, plot(regul, 'LineWidth',1.5)
title('Regularization term')
axis tight
title(t,['Iterations, \lambda=',num2str(lambda,2),...
    ', \tau=',num2str(tau,2)])

%% Selecting ROI
swsPath = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\TV\sws';
load([swsPath,'\1.mat']);
x0inc = 15; z0 = 11; L = 11; x0back = 1.5;
figure, 
imagesc(dinf.x*1000,dinf.z*1000, Bmode); colormap gray
colorbar
axis equal
hold on
rectangle('Position',[x0inc z0 L L], 'LineWidth',2),
rectangle('Position',[x0back z0 L L], 'LineWidth',2, 'EdgeColor','w'),
hold off
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title('B-mode ROI')
ax = gca; ax.FontSize = 12;

%% Calculating CNR
[X,Z] = meshgrid(dinf.x*1000,dinf.z*1000);
maskInc = (X>x0inc & X<x0inc+L & Z>z0 & Z<z0+L);
maskBack = (X>x0back & X<x0back+L & Z>z0 & Z<z0+L);

cnrTV = NaN(1,Nim);
cnrRW = NaN(1,Nim);
cnrTNV = NaN(1,Nim);
cnrCWT = NaN(1,Nim);

meanIncTV = NaN(1,Nim);
meanBackTV = NaN(1,Nim);
stdIncTV = NaN(1,Nim);
stdBackTV = NaN(1,Nim);

meanIncRW = NaN(1,Nim);
meanBackRW = NaN(1,Nim);
stdIncRW = NaN(1,Nim);
stdBackRW = NaN(1,Nim);

meanIncTNV = NaN(1,Nim);
meanBackTNV = NaN(1,Nim);
stdIncTNV = NaN(1,Nim);
stdBackTNV = NaN(1,Nim);

meanIncCWT = zeros(1,Nim);
meanBackCWT = zeros(1,Nim);
stdIncCWT = zeros(1,Nim);
stdBackCWT = zeros(1,Nim);
for iCuac = 1:Nim
    iIm = selectedImages(iCuac);
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

    %swsTNVchannel = swsTNV(:,:,iIm);
    swsTNVchannel = swsTNV(:,:,iCuac);
    swsInc = swsTNVchannel(maskInc);
    swsBack = swsTNVchannel(maskBack);
    cnrTNV(iIm) = 2*(mean(swsBack) - mean(swsInc))^2 / ...
        (var(swsInc) + var(swsBack));
    meanIncTNV(iIm) = mean(swsInc);
    meanBackTNV(iIm) = mean(swsBack);
    stdIncTNV(iIm) = std(swsInc);
    stdBackTNV(iIm) = std(swsBack); 

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
VibFreqArray = 200:50:450;
figure('Position', [100 100 600 400]),
plot(VibFreqArray,db(cnrRW),'o-', 'LineWidth',2)
hold on
plot(VibFreqArray,db(cnrTV),'o-', 'LineWidth',2)
plot(VibFreqArray,db(cnrTNV),'o-', 'LineWidth',2)
hold off
legend({'R-WAVE','TV','TNV'}, 'Location','northwest');
ylabel('CNR [dB]'), xlabel('F_v')
grid on
xlim([180 470])

%% Comparing mean SWS and std
figure('Position', [100 100 400 400]),
tiledlayout(1,4)
nexttile
errorbar(VibFreqArray,meanBackCWT,stdBackCWT, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncCWT,stdIncCWT, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 470]), ylim([2 4])
title('CWT')

nexttile
errorbar(VibFreqArray,meanBackRW,stdBackRW, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncRW,stdIncRW, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 470]), ylim([2 4])
title('R-WAVE')

nexttile
errorbar(VibFreqArray,meanBackTV,stdBackTV, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncTV,stdIncTV, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 470]), ylim([2 4])
title('Total Variation')

nexttile
errorbar(VibFreqArray,meanBackTNV,stdBackTNV, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncTNV,stdIncTNV, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 470]), ylim([2 4])
title('TNV')

%% Comparing lateral profile
iz = 70;
freqArr = 200:20:360;
for iIm = 1:Nim
    load([swsPath,'\',num2str(iIm),'.mat']);
    x = dinf.x*1000;
    lineRW = swsRW(iz,:);
    lineTV = swsTV(iz,:);
    lineTNV = swsTNV(iz,:,iIm);

    figure('Position', [100 100 600 250]),
    plot(x,lineRW, 'LineWidth',1.5)
    hold on
    plot(x,lineTV, 'LineWidth',1.5)
    plot(x,lineTNV, 'LineWidth',1.5)
    hold off
    legend({'RW','TV','TNV'})
    axis tight
    ylim([3 5.5])
    zi = dinf.z(iz)*1000;
    title(['Lateral profile at z=',num2str(zi,2),'mm, f_v=',...
        num2str(freqArr(iIm)),'Hz'])
end
%% Comparing axial profile
ix = 80;
freqArr = 200:20:360;
for iIm = 1:Nim
    load([swsPath,'\',num2str(iIm),'.mat']);
    z = dinf.z*1000;
    lineRW = swsRW(:,ix);
    lineTV = swsTV(:,ix);
    lineTNV = swsTNV(:,ix,iIm);

    figure('Position', [100 100 600 250]),
    plot(z,lineRW, 'LineWidth',1.5)
    hold on
    plot(z,lineTV, 'LineWidth',1.5)
    plot(z,lineTNV, 'LineWidth',1.5)
    hold off
    legend({'RW','TV','TNV'})
    axis tight
    ylim([3 6])
    xi = dinf.x(ix)*1000;
    title(['Axial profile at x=',num2str(xi,2),'mm, f_v=',...
        num2str(freqArr(iIm)),'Hz'])
end

%% Comparing lateral resolution
iz = 70;
ix0 = 23; ixf = 80;
%ix0 = 23; ixf = 110;
ix0n = 70; ixfn = 110;
freqArr = 200:20:360;
fwhmRW = NaN(1,Nim);
fwhmTV = NaN(1,Nim);
fwhmTNV = NaN(1,Nim);

figure('Position', [100 100 1200 600]),
tiledlayout(3,3)
for iCuac = 1:Nim
    iIm = selectedImages(iCuac);
    load([swsPath,'\',num2str(iIm),'.mat']);
    x = dinf.x(ix0:ixf-1)*1000;
    lineRW = diff(swsRW(iz,ix0:ixf));
    lineTV = diff(swsTV(iz,ix0:ixf));
    lineTNV = diff(swsTNV(iz,ix0:ixf,iCuac));
    lineRWneg = diff(swsRW(iz,ix0n:ixfn));
    lineTVneg = diff(swsTV(iz,ix0n:ixfn));
    lineTNVneg = diff(swsTNV(iz,ix0n:ixfn,iCuac));
    fwhmRW(iIm) = fwhm(x,lineRW) + fwhm(x,-lineRWneg);
    fwhmTV(iIm) = fwhm(x,lineTV) + fwhm(x,-lineTVneg);
    fwhmTNV(iIm) = fwhm(x,lineTNV) + fwhm(x,-lineTNVneg);
    nexttile,
    plot(lineRW, 'LineWidth',1.5)
    hold on
    plot(lineTV, 'LineWidth',1.5)
    plot(lineTNV, 'LineWidth',1.5)
    hold off
    legend({'RW','TV','TNV'}, 'NumColumns',3, 'Location','north')
    axis tight
    ylim([-.4 .4])
    zi = dinf.z(iz)*1000;
    title(['Lateral difference at z=',num2str(zi,2),'mm, f_v=',...
        num2str(freqArr(iIm)),'Hz'])
end
%%
figure('Position', [100 100 600 400]),
plot(VibFreqArray,fwhmRW,'o-', 'LineWidth',2)
hold on
plot(VibFreqArray,fwhmTV,'o-', 'LineWidth',2)
plot(VibFreqArray,fwhmTNV,'o-', 'LineWidth',2)
hold off
legend({'R-WAVE','TV','TNV'}, 'Location','northwest');
ylabel('FWHM [mm]'), xlabel('F_v')
grid on
xlim([180 470])
ylim([2 14])

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