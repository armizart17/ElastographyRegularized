% Script para procesar data de adquisision WIDE BEAM 20-09
%% Loading filtered sono signal
clear, clc
addpath('./RegularizationFunctions/',"./ElastographyFunctions")

baseDir = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\20_09_23';
rawDir = [baseDir,'\raw'];
rawFiles = dir([rawDir,'\*.mat']);
Nim = length(rawFiles);

iqDir = [baseDir,'\iq'];
sonoDir = [baseDir,'\sono'];
swsDir = [baseDir,'\sws'];

% Colormap
load('MyColormaps.mat')

%% Processing raw files
if ~exist(iqDir,'dir'), mkdir(iqDir); end

for iAcq = 1:Nim
    % iAcq = 6;
        file = rawFiles(iAcq);
        load([file.folder,'\',file.name]);
    
    % Acquisition parameters
    framePeriod = (P.numRays*SeqControl(2).argument + SeqControl(3).argument + ...
        SeqControl(6).argument*P.PRIsDop*numRaysDop + SeqControl(8).argument)*1e-6;
    dinf.frameRate = 1/framePeriod; % ta bien, un ciclo son 30 frames
    dinf.vibf = str2double(file.name(5:7));
    dinf.offset = 2; % MANUALLY SET
    dinf.PRF = P.PRFDop;
    dinf.fc = Trans.frequency*1e6;
    dinf.c0 = Resource.Parameters.speedOfSound;
    dinf.wl = dinf.c0/dinf.fc;
    dinf.fs = 2*dinf.fc; % No se
    
    % Bmode data
    IQ = IData{1} + 1j*QData{1};
    dinf.Bmode = db(IQ);
    
    % Bmode pixel data
    dinf.dxBmode = PData(1).PDelta(1)*dinf.wl;
    dinf.dzBmode = PData(1).PDelta(3)*dinf.wl;
    x0 = PData(1).Origin(1)*dinf.wl;
    z0 = PData(1).Origin(3)*dinf.wl;
    dinf.xBmode = (0:PData(1).Size(2)-1)*dinf.dxBmode + x0;
    dinf.zBmode = (0:PData(1).Size(1)-1)*dinf.dzBmode + z0;
    
    % Doppler IQ Data
    IQ = squeeze(IData{2} + 1j*QData{2});
    [Nz,Nx,Nens,Nframes] = size(IQ);
    
    % Doppler pixel data
    dinf.dx = PData(2).PDelta(1)*dinf.wl;
    dinf.dz = PData(2).PDelta(3)*dinf.wl;
    x0 = PData(2).Origin(1)*dinf.wl;
    z0 = PData(2).Origin(3)*dinf.wl;
    dinf.x = (0:Nx-1)*dinf.dx + x0;
    dinf.z = (0:Nz-1)*dinf.dz + z0;
    
    % Cropping ROI
    mask = zeros(Nz,Nx);
    mask(PData(2).Region(5).PixelsLA) = 1;
    %figure, imagesc(dinf.x*1000,dinf.z*1000,mask)
    dinf.x = dinf.x(logical(mask(50,:)));
    dinf.z = dinf.z(logical(mask(:,50)));
    maskIQ = logical(repmat(mask,[1,1,Nens,Nframes]));
    IQ = reshape(IQ(maskIQ),length(dinf.z),length(dinf.x),Nens,Nframes);
    IQ = IQ(1:end-1,:,:,:);
    dinf.z = dinf.z(1:end-1);
    
    % Plotting B-mode images
    % figure, tiledlayout(1,2);
    % nexttile,
    % x = dinf.xBmode*1000;
    % z = dinf.zBmode*1000;
    % imagesc(x,z,dinf.Bmode)
    % title('Bmode')
    % axis equal
    % axis tight
    % colormap gray
    % 
    % nexttile
    % x = dinf.x*1000;
    % z = dinf.z*1000;
    % imagesc(x,z,db(IQ(:,:,1,1)))
    % title('Doppler area')
    % axis equal
    % axis tight
    % colormap gray
    save([iqDir,'\',num2str(iAcq),'.mat'],'IQ','dinf')
    
    % Generating sono
    sono = generateSonoFromIQ(IQ,dinf.PRF);   % Sonoelastography
    
    Prop.FrameRate = dinf.frameRate; Prop.pitch = dinf.dx;
    Prop.VibFreq = dinf.vibf; Prop.VibFreqOffset = 2;
    move = 'right'; swsRange = [2,7];
    [~,sonoFilt,~] = process_sono_data(sono,Prop,move,swsRange);
    
    % Plotting to obtain start index
    signal = squeeze(sonoFilt(50,50,:));
    figure, plot(signal)
    title(['Sono line F_v=',num2str(dinf.vibf),'Hz'])
    axis tight
end
%% Generating aligned sono
close all,
if ~exist(sonoDir,'dir'), mkdir(sonoDir); end
for iIm = 1:Nim
    load([iqDir,'\',num2str(iIm),'.mat'])
    % Generating sono
    sono = generateSonoFromIQ(IQ,dinf.PRF);   % Sonoelastography
    
    Prop.FrameRate = dinf.frameRate; Prop.pitch = dinf.dx;
    Prop.VibFreq = dinf.vibf; Prop.VibFreqOffset = 2;
    move = 'right'; swsRange = [2,7];
    
    iStart = [24,35,40,20,44,16]; % OBTENIDOS MANUALMENTE
    sonoNew = cat(3,sono(:,:,iStart(iIm):end),sono(:,:,1:iStart(iIm)-1));
    [sonoFiltMov,sonoFilt,~] = process_sono_data(sonoNew,Prop,move,swsRange);
    Nframes = size(sono,3);
    figure,
    for iFrame = 1:Nframes
        x = dinf.x*1000;
        z = dinf.z*1000;
        imagesc(x,z,sonoFiltMov(:,:,iFrame), 1*[-1 1])
        axis equal
        axis tight
        colormap(sonomap)
        title(['Filtered sono F_v=',num2str(dinf.vibf),'Hz'])
        pause(1/dinf.frameRate)
    end
    sono = sonoFiltMov;
    save([sonoDir,'\',num2str(iIm),'.mat'],'sono','dinf')
end

%% Plotting Bmode and sono
fig1 = figure('Position',[100 100 1000 500]);
t1 = tiledlayout(fig1,2,3);
fig2 = figure('Position',[100 100 1000 500]);
t2 = tiledlayout(fig2,2,3);
for iIm = 1:Nim
    load([sonoDir,'\',num2str(iIm),'.mat'])

    x = dinf.x*1000;
    z = dinf.z*1000;
    nexttile(t1,iIm);
    imagesc(x,z,sono(:,:,1),1.5*[-1 1])
    colormap(sonomap)
    colorbar
    axis equal
    axis tight
    title(['Sono F_v=',num2str(dinf.vibf)])
    
    xBm = dinf.xBmode*1000;
    zBm = dinf.zBmode*1000;
    nexttile(t2,iIm);
    imagesc(xBm,zBm,dinf.Bmode, [110 170])
    colormap gray
    colorbar
    axis equal
    xlim([x(1),x(end)]), ylim([z(1),z(end)])
    title(['B-mode F_v=',num2str(dinf.vibf)])
end

%% Independent SWS maps
close all,

if ~exist(swsDir,'dir'), mkdir(swsDir); end

ParamsTV.mu = 8;
ParamsTV.tol = 1e-1;
ParamsTV.isotropic = true;

RW_Param.dual = boolean(1); RW_Param.k = 1;
RW_Param.N_window = 20; RW_Param.beta = 1/100000;
RW_Param.tolerance = 1/1000;RW_Param.operator = 'G';
RW_Param.alpha = 2;

fig1 = figure('Position',[100 100 1000 500]);
t1 = tiledlayout(fig1,2,3);
fig2 = figure('Position',[100 100 1000 500]);
t2 = tiledlayout(fig2,2,3);
fig3 = figure('Position',[100 100 1000 500]);
t3 = tiledlayout(fig3,2,3);
for iIm = 1:Nim
    load([sonoDir,'\',num2str(iIm),'.mat'])

    % fprintf("\nVibration Frequency = %d Hz\n",dinf.vibf);
    % Properties.pitch = dinf.dx;
    % Properties.VibFreq = dinf.vibf;
    % [swsTV,C] = calculateSWSTV(sono,Properties,ParamsTV);
    % fprintf("Number of iterations: %d\n",length(C));
    % % 
    % tic
    % swsRW = rWave2(sono,Properties,RW_Param);
    % t = toc;
    % fprintf('Exec. time for R-W: %f\n',t)
    % 
    % tic
    % swsCWT = process_CWT(sono,Properties,[1 16]);
    % t = toc;
    % fprintf('Exec. time for CWT: %f\n',t)
    % 
    % save([swsDir,'\',num2str(iIm),'.mat'],...
    %     'swsTV','C','swsRW','swsCWT','dinf');

    % Selecting ROI
    x = dinf.x*1000;
    z = dinf.z*1000;
    SWS_im_range = [2 6];
    
    nexttile(t1,iIm);
    im = imagesc(x,z,swsTV,SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['SWS from TV, \mu=',num2str(ParamsTV.mu,2),...
        ', f_v=',num2str(dinf.vibf)])
    ax = gca; ax.FontSize = 12;

    nexttile(t2,iIm);
    im = imagesc(x,z,swsRW,SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['SWS from RW, \alpha=',num2str(RW_Param.alpha,2),...
        ', f_v=',num2str(dinf.vibf)])
    ax = gca; ax.FontSize = 12;
 
    nexttile(t3,iIm);
    swsCwtIm = squeeze(mean(swsCWT,3));
    im = imagesc(x,z,swsCwtIm,SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['SWS from CWT, f_v=',num2str(dinf.vibf)])
    ax = gca; ax.FontSize = 12;
end

%% Selecting ROI
x0inc = -2; z0 = 14; L = 9; x0back = -15;
x = dinf.x*1000;
z = dinf.z*1000;
figure, 
imagesc(dinf.xBmode*1000,dinf.zBmode*1000, dinf.Bmode); colormap gray
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
    load([swsDir,'\',num2str(iIm),'.mat']);

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
figure('Position', [100 100 1000 400]),
tiledlayout(1,3)
nexttile
errorbar(VibFreqArray,meanBackTV,stdBackTV, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncTV,stdIncTV, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 470]), ylim([2 6])
title('Total Variation')

nexttile
errorbar(VibFreqArray,meanBackRW,stdBackRW, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncRW,stdIncRW, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 470]), ylim([2 6])
title('R-WAVE')

nexttile
errorbar(VibFreqArray,meanBackCWT,stdBackCWT, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncCWT,stdIncCWT, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 470]), ylim([2 6])
title('CWT')

%% --------------------------------------------------------------------- %%
%% TOTAL NUCLEAR VARIATION
%selectedImages = [3,4,5,7,8,9];
selectedImages = 1:6;
Nim = length(selectedImages);
for iIm = 1:Nim
    % Loading
    load([sonoDir,'\',num2str(selectedImages(iIm)),'.mat'])
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
for lambda = logspace(-0.5,1.5,5)
    for tau = logspace(-1.5,0.5,5)
        %%
        % lambda = 1; % +-50%
        % tau = 0.2;   % +-50%
        maxIter = 500;
        tol = 3e-4;
        stableIter = 50;
        Nim = length(selectedImages);
        tic
        [u, cost, error, fide, regul] = pdo_inv_tnv(B, Nz, Nx, A, ...
            lambda, tau, maxIter, tol, Nim, stableIter);
        swsTNV = reshape(u,Nz,Nx,Nim);
        toc
        %save([swsDir,'\TNV.mat'],'swsTNV','dinf');
        
        % Results
        %load([swsDir,'\TNV.mat'],'swsTNV');
        x = dinf.x * 1000;
        z = dinf.z * 1000;
        SWS_im_range = [2,6];
        VibFreqArray = 200:20:360;
        
        figure('Position',[100 100 1000 500]);
        %t = tiledlayout(3,3);
        t = tiledlayout(2,3);
        for iIm = 1:Nim
            nexttile
            im = imagesc(x,z,squeeze(swsTNV(:,:,iIm)),SWS_im_range);
            %im.AlphaData = dinf.valid;
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
    end
end
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




