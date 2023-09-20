% Script para procesar data del drive de Arroyo (mama y phantoms)
%% Loading filtered sono signal
clear, clc
addpath('./RegularizationFunctions/',"./ElastographyFunctions")

baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets' ...
    '\VibroElastografia\Inclusion\06-23-2017-Generic'];

rawDir = [baseDir,'\RF'];
crfFiles = dir([rawDir,'\*.crf']);
rfFiles = dir([rawDir,'\*.rf']);
sonoDir = [baseDir,'\sono'];
croppedDir = [baseDir,'\cropped'];
swsDir = [baseDir,'\sws'];

% Colormap
load('MyColormaps.mat')

%% Getting sono signals
if ~exist(sonoDir,'dir'), mkdir(sonoDir); end
wf = [1,0];
Nim = length(crfFiles);
for iFile = 1:Nim
    name_file = crfFiles(iFile).name;                % Name of CRF file
    name_file2 = rfFiles(iFile).name;               % Name of RF file
    f = str2double(name_file(1:3));          % Vibration frequency
    
    file = [rawDir,'\',name_file];        % Path of the CRF file
    file2 = [rawDir,'\',name_file2];      % Path of the RF file
    [Data, Properties_CRF] = RPread(file);     % Acquiring CRF data
    [Data2,Properties_RF] = RPread(file2);   % Acquiring RF data

    if isempty(Data) || isempty(Data2)
        continue
    else
        Properties.PRF = Properties_CRF.dr;          % PRF
        Properties.FrameRate = Properties_RF.dr;  % Frame rate
        IQ = process_IQ(Data,Properties_CRF);        % IQ demodulation of CRF
        sono = process_colorflowIQ_m(IQ,Properties_CRF,wf);   % Sonoelastography
        
        Properties.pitch = 3.08e-04;
        Depth = 1540/(2*Properties_CRF.sf) * Properties_CRF.h;
        Properties.z = linspace(0,Depth,Properties_CRF.h/5);
        Properties.x = (1:Properties_CRF.w) * Properties.pitch;
        
        Properties.VibFreq = f;
        Properties.VibFreqOffset = 0.4;

        Bmode = db(squeeze(IQ(:,1,:,1))); 
        Bmode = Bmode - max(Bmode(:));
        Properties.Bmode = Bmode;
        
        save([sonoDir,'\', num2str(iFile)],'sono','Properties','Bmode')
    end
end

%% pre-filtering sono
swsRange = [2,8];  
move = 'left';
[~,sonoFilt,~] = process_sono_data(sono,Properties,move,swsRange);
signal = squeeze(sonoFilt(50,50,:));
Nframes = size(sono,3);
t = (0:Nframes-1)/Properties.FrameRate;
f = (0:Nframes-1)/Nframes * Properties.FrameRate;
figure, tiledlayout(2,1)
nexttile
plot(t,signal)
nexttile
plot(f,abs(fft(signal)))

%% Subsampling and cropping
if ~exist(croppedDir,'dir'), mkdir(croppedDir); end

% Filtering
swsRange = [2,8];  
move = 'left';

% Cropping
z0 = 0e-3;
zf = 25e-3;

% Looping
fig1 = figure('Position',[100 100 1000 500]);
t1 = tiledlayout(fig1,3,3);
fig2 = figure('Position',[100 100 1000 500]);
t2 = tiledlayout(fig2,3,3);
for iIm = 1:Nim
    % Pre-processing
    load([sonoDir,'\', num2str(iIm)]);
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
    zNew = Properties.z(1:2:end);
    izROI = zNew>z0 & zNew<zf;
    sono = sonoNew(izROI,:,:);
    Properties.z = zNew(izROI);
    
    x = Properties.x*1000;
    z = Properties.z*1000;
    nexttile(t1,iIm);
    imagesc(x,z,sono(:,:,1),1.5*[-1 1])
    colormap(sonomap)
    colorbar
    axis equal
    axis tight
    title(['Sono f_v=',num2str(Properties.VibFreq)])

    nexttile(t2,iIm);
    imagesc(x,z,Properties.Bmode)
    colormap gray
    colorbar
    axis equal
    axis tight
    title(['B-mode f_v=',num2str(Properties.VibFreq)])
    
    save([croppedDir,'\',num2str(iIm),'.mat'],"sono","Properties");
end

%% --------------------------------------------------------------------- %%
%% 2-DIMENSIONAL TOTAL VARIATION
Nim = 8; % Number of channels
if ~exist(swsDir,'dir'), mkdir(swsDir); end
ParamsTV.mu = 40;
ParamsTV.tol = 1e-1;
ParamsTV.isotropic = true;

RW_Param.dual = boolean(1); RW_Param.k = 1;
RW_Param.N_window = 20; RW_Param.beta = 1/100000;
RW_Param.tolerance = 1/1000;RW_Param.operator = 'G';
RW_Param.alpha = 5;

fig1 = figure('Position',[100 100 1000 500]);
t1 = tiledlayout(fig1,3,3);
fig2 = figure('Position',[100 100 1000 500]);
t2 = tiledlayout(fig2,3,3);
fig3 = figure('Position',[100 100 1000 500]);
t3 = tiledlayout(fig3,3,3);

for iIm = 1:Nim
    load([croppedDir,'\',num2str(iIm),'.mat'])
    load([swsDir,'\',num2str(iIm),'.mat']);

%     fprintf("\nVibration Frequency = %d Hz\n",Properties.VibFreq);
%     [swsTV,C] = calculateSWSTV(sono,Properties,ParamsTV);
%     fprintf("Number of iterations: %d\n",length(C));
%     % 
    tic
    swsRW = rWave2(sono,Properties,RW_Param);
    t = toc;
    fprintf('Exec. time for R-W: %f\n',t)
    
%     tic
%     swsCWT = process_CWT(sono,Properties,[1 16]);
%     t = toc;
%     fprintf('Exec. time for CWT: %f\n',t)
    
    save([swsDir,'\',num2str(iIm),'.mat'],...
        'swsTV','C','swsRW','swsCWT','Properties');

    % Selecting ROI
    x = Properties.x*1000;
    z = Properties.z*1000;
    SWS_im_range = [4 7];
    
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
Nim = 8;
for iIm = 1:Nim
    % Loading
    load([croppedDir,'\',num2str(iIm),'.mat'])
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
% larger tau increases num of iterations and makes the results less
% dependendent on lambda
% smaller tau increases agreement between channels
lambda = 40; % 5 +-50%
tau = 0.2;   % 0.5 +-50%
maxIter = 1000;
tol = 3e-4;
stableIter = 50;
[u, cost, error, fide, regul] = pdo_inv_tnv(B, Nz, Nx, A, ...
    lambda, tau, maxIter, tol, Nim, stableIter);
swsTNV = reshape(u,Nz,Nx,Nim);
save([swsDir,'\TNV.mat'],'swsTNV');

%% Results
load([swsDir,'\TNV.mat'],'swsTNV');
x = Properties.x * 1000;
z = Properties.z * 1000;
SWS_im_range = [4,7];
VibFreqArray = 200:20:360;

figure('Position',[100 100 1200 600]);
t = tiledlayout(3,3);
for iIm = 1:Nim
    nexttile
    imagesc(x,z,squeeze(swsTNV(:,:,iIm)),SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['F_v = ',num2str(VibFreqArray(iIm)),' Hz'])
    ax = gca; ax.FontSize = 12;
end
title(t,['SWS from TNV, \lambda=',num2str(lambda,2),...
    ', \tau=',num2str(tau,2)])

figure,
t = tiledlayout(2,2);
nexttile, plot(cost, 'LineWidth',1.5)
title('Cost function')
xlim([0 120])
nexttile, plot(error, 'LineWidth',1.5)
title('Relative error')
xlim([0 120])
nexttile, plot(fide, 'LineWidth',1.5)
title('Fidelity term')
xlim([0 120])
nexttile, plot(regul, 'LineWidth',1.5)
title('Regularization term')
xlim([0 120])
title(t,['Iterations, \lambda=',num2str(lambda,2),...
    ', \tau=',num2str(tau,2)])

%% Selecting ROI
load([swsDir,'\1.mat']);
x0inc = 15; z0 = 11; L = 11; x0back = 1.5;
figure, 
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
title('B-mode ROI')
ax = gca; ax.FontSize = 12;

%% Calculating CNR
VibFreqArray = 200:20:360;
[X,Z] = meshgrid(Properties.x*1000,Properties.z*1000);
maskInc = (X>x0inc & X<x0inc+L & Z>z0 & Z<z0+L);
maskBack = (X>x0back & X<x0back+L & Z>z0 & Z<z0+L);

cnrTV = zeros(1,Nim);
cnrRW = zeros(1,Nim);
cnrTNV = zeros(1,Nim);

meanIncTV = zeros(1,Nim);
meanBackTV = zeros(1,Nim);
stdIncTV = zeros(1,Nim);
stdBackTV = zeros(1,Nim);

meanIncRW = zeros(1,Nim);
meanBackRW = zeros(1,Nim);
stdIncRW = zeros(1,Nim);
stdBackRW = zeros(1,Nim);

meanIncTNV = zeros(1,Nim);
meanBackTNV = zeros(1,Nim);
stdIncTNV = zeros(1,Nim);
stdBackTNV = zeros(1,Nim);
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

    swsTNVchannel = swsTNV(:,:,iIm);
    swsInc = swsTNVchannel(maskInc);
    swsBack = swsTNVchannel(maskBack);
    cnrTNV(iIm) = 2*(mean(swsBack) - mean(swsInc))^2 / ...
        (var(swsInc) + var(swsBack));
    meanIncTNV(iIm) = mean(swsInc);
    meanBackTNV(iIm) = mean(swsBack);
    stdIncTNV(iIm) = std(swsInc);
    stdBackTNV(iIm) = std(swsBack); 
end
%% Comparing CNR
figure('Position', [100 100 600 400]),
plot(VibFreqArray,db(cnrRW),'o-', 'LineWidth',2)
hold on
plot(VibFreqArray,db(cnrTV),'o-', 'LineWidth',2)
plot(VibFreqArray,db(cnrTNV),'o-', 'LineWidth',2)
hold off
legend({'R-WAVE','TV','TNV'}, 'Location','northwest');
ylabel('CNR [dB]'), xlabel('F_v')
grid on
xlim([180 380])

%% Comparing mean SWS and std
figure('Position', [100 100 400 400]),
tiledlayout(1,3)
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
errorbar(VibFreqArray,meanBackTNV,stdBackTNV, 'LineWidth',2)
hold on
errorbar(VibFreqArray,meanIncTNV,stdIncTNV, 'LineWidth',2)
hold off
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('Vibration frequency [Hz]')
grid on
xlim([180 380])
title('TNV')

%% Comparing lateral profile
iz = 70;
freqArr = 200:20:360;
for iIm = 1:Nim
    load([swsPath,'\',num2str(iIm),'.mat']);
    x = Properties.x*1000;
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
    zi = Properties.z(iz)*1000;
    title(['Lateral profile at z=',num2str(zi,2),'mm, f_v=',...
        num2str(freqArr(iIm)),'Hz'])
end
%% Comparing axial profile
ix = 80;
freqArr = 200:20:360;
for iIm = 1:Nim
    load([swsPath,'\',num2str(iIm),'.mat']);
    z = Properties.z*1000;
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
    xi = Properties.x(ix)*1000;
    title(['Axial profile at x=',num2str(xi,2),'mm, f_v=',...
        num2str(freqArr(iIm)),'Hz'])
end

%% Comparing lateral resolution
iz = 70;
ix0 = 23; ixf = 80;
%ix0 = 23; ixf = 110;
ix0n = 70; ixfn = 110;
freqArr = 200:20:360;
fwhmRW = zeros(1,Nim);
fwhmTV = zeros(1,Nim);
fwhmTNV = zeros(1,Nim);

figure('Position', [100 100 1200 600]),
tiledlayout(3,3)
for iIm = 1:Nim
    load([swsPath,'\',num2str(iIm),'.mat']);
    x = Properties.x(ix0:ixf-1)*1000;
    lineRW = diff(swsRW(iz,ix0:ixf));
    lineTV = diff(swsTV(iz,ix0:ixf));
    lineTNV = diff(swsTNV(iz,ix0:ixf,iIm));
    lineRWneg = diff(swsRW(iz,ix0n:ixfn));
    lineTVneg = diff(swsTV(iz,ix0n:ixfn));
    lineTNVneg = diff(swsTNV(iz,ix0n:ixfn,iIm));
    fwhmRW(iIm) = fwhm(x,lineRW) + fwhm(x,-lineRWneg);
    fwhmTV(iIm) = fwhm(x,lineTV) + fwhm(x,-lineTVneg);
    fwhmTNV(iIm) = fwhm(x,lineTNV) + fwhm(x,-lineTNVneg);
    nexttile,
    plot(lineRW, 'LineWidth',1.5)
    hold on
    plot(lineTV, 'LineWidth',1.5)
    plot(lineTNV, 'LineWidth',1.5)
    hold off
    legend({'RW','TV','TNV'}, 'NumColumns',3)
    axis tight
    ylim([-.4 .4])
    zi = Properties.z(iz)*1000;
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
xlim([180 380])

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