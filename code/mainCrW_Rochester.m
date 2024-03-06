%% THIS MODIFICATION WAS MADE BY EMZ (March 2024)
% Test for new data Juvenal Ormachea, PhD.
% University of Rochester

% Hello Prof. Parker
% 
% I could find some CrW experiment files. Here are 60 different experiments
% that I collected in a cone phantom. Twenty different vibration frequencies
% at three different positions.
% 
% I think the only thing missing are the dimension parameters. 
% Roberto's student can open each correspondent .mat file and run these 
% lines for the axial and lateral sizes in mm, respectively:
% sa 
% zz=10^3*d1.movies.simplexMovies(1,3).params(128,1).value;
% xx=10^3*d1.movies.simplexMovies(1,3).params(133,1).value;
% 
% For the frame rate it should be:
% FR=d1.movies.simplexMovies(1,3).params(137,1).value;
% 
% I hope these files are good for the PUCP people.
% All the best
% Juve
% 
% PS. Because the big size of this folder, I'll keep it in my share drive 
% just for a week. Please, ask them to confirm with you when they 
% have download the data.
% 

% NOTES EMZ March
% Image 1-20  120Hz - 500Hz
% Image 21-40 120Hz - 500Hz
% Image 41-60 120Hz - 500Hz

%% Loading filtered sono signal
clear, clc

% addpath('./data/ForRoberto/')
addpath('./RegularizationFunctions/',"./ElastographyFunctions/", "./utils/")

baseDir = './'; % CHANGE ACCORDING TO WHAT YOU WANT
dataDir = [baseDir,'data/ForRoberto/09-26-14/'];


resultsDir =  './out/ForRoberto/';
sonoPath = [resultsDir,'sono/'];
swsPath = [resultsDir,'sws/'];
outiniPath = [resultsDir,'iniResults/'];

% if ~exist("resultsDir","dir"); mkdir(resultsDir); end
% if ~exist("sonoPath","dir"); mkdir(sonoPath); end
% if ~exist("swsPath","dir"); mkdir(swsPath); end
% if ~exist("outiniPath","dir"); mkdir(outiniPath); end

% Colormap
load('./MyColormaps.mat')
%%
numCase = 1;
% OPTION 1
% st_sono = load([dataDir, 'MatlabProcessed/Image', num2str(numCase),'/sono.mat']);
% st_sono_filt = load([dataDir, 'MatlabProcessed/Image', num2str(numCase),'/sono_filt.mat']);
% OPTION 2
% load(dataDir + "MatlabProcessed/Image" + num2str(numCase) + "/sono.mat")
% load(dataDir + "MatlabProcessed/Image" + num2str(numCase) + "/sono_filt.mat")


load([dataDir, 'MatlabProcessed/Image', num2str(numCase),'/sono.mat']);
load([dataDir, 'MatlabProcessed/Image', num2str(numCase),'/sono_filt.mat']);

%%% my modification 2

% st_E9LJP = load([dataDir, 'E9LJPNO0_DCMSTRUCT.mat']);
depth=1e3*d1.movies.simplexMovies(1,3).params(128,1).value;
width=1e3*d1.movies.simplexMovies(1,3).params(133,1).value;


%%


fig1 = figure('Position',[100 100 1000 500]);
t1 = tiledlayout(fig1,1,3);

figure, 
subplot(1,3,1), 
imagesc(x_b, z_b, Bmode), title('Bmode')
xlabel('Lateral [mm]'), ylabel('Axial[mm]')

subplot(1,3,2), 
imagesc(x_s, z_s, sono(:,:,1)), title('Sono')
xlabel('Lateral [mm]'), ylabel('Axial[mm]')


subplot(1,3,3), 
imagesc(x_s, z_s, sono_filt(:,:,1)), title('Sono Filt')
xlabel('Lateral [mm]'), ylabel('Axial[mm]')
%% Subsampling and cropping
% Filtering
swsRange = [2,8]; % [m/s]
move = 'left';


%%%%%%%%%%%%%%%%% AXIS extraction Attempt 1 %%%%%%%%%%%%%%%%%
numCase = 1;
load([dataDir, 'MatlabProcessed/Image', num2str(numCase),'/sono.mat']);
load([dataDir, 'MatlabProcessed/Image', num2str(numCase),'/sono_filt.mat']);

depth = 60e-3;     % [m]
width = 38.4e-3;   % [m]
[M, N, P] = size(sono); 

% B mode axis
z_b = linspace(0, depth, size(Bmode,1));
x_b = linspace(-width/2, width/2, size(Bmode,2));

% Sono axis
z_s = linspace(0, depth, M);
x_s = linspace(-width/2, width/2, N);

%%%%%%%%%%%%%%%%% AXIS extraction Attempt 1 %%%%%%%%%%%%%%%%%

% Cropping
z0 = 5e-3;  % [m]
zf = 50e-3; % [m]

% Looping
% Nacq = input('Enter # of acquisitions: ');
Nacq = 3; % Number of acquisitions
% Nim = input('Enter # of freq. channels: ');
Nim = 20; % Number of channels
VibFreqArray = 120:20:500; % [Hz]

% Plots
font = 14;
caxis_bmode = [-60 0]; % [dB]


for iacq = 1:Nacq

fig1 = figure(10*iacq);
set(10*iacq, 'Position',[100 100 1500 1000]);
t1 = tiledlayout(fig1,4,5);
sgtitle(['\bf Sono Filt Acq #', num2str(iacq)], 'FontSize', font );

fig2 = figure(10*iacq+1);
set(10*iacq+1, 'Position',[100 100 1500 1000]);
t2 = tiledlayout(fig2,4,5);
sgtitle(['\bf B-mode Acq #', num2str(iacq)], 'FontSize', font );

for iIm = 1:Nim
    % Pre-processing

    load([dataDir, 'MatlabProcessed/Image', num2str(10*iacq+iIm),'/sono.mat']);
    load([dataDir, 'MatlabProcessed/Image', num2str(10*iacq+iIm),'/sono_filt.mat']);


    % Subsampling
    R = 2; % DECIMATION FACTOR
    [Nz,Nx,Nt] = size(sono_filt);
    sonoNew = zeros([ceil(Nz/R),Nx,Nt]);
    [b,a] = butter(4,1/R);
    for ix = 1:Nx
        for it = 1:Nt
            signal = filtfilt(b,a,sono_filt(:,ix,it));
            sonoNew(:,ix,it) = signal(1:R:end);
        end
    end
    
    % Selecting ROI
    zNew = z_s(1:2:end);
    izROI = zNew>z0 & zNew<zf;
    sono = sonoNew(izROI,:,:);
    Properties.Depth_S = zNew(izROI);
    izBmode = z_b >z0 & z_b <zf;
    Properties.Bmode = Bmode(izBmode,:);
    Properties.Depth_B = z_b(izBmode);
    Properties.Width_B = x_b;
    
    %%%%%%%%%%%%%% PLOT PATTERN %%%%%%%%%%%%%%
    x = x_s*1e3; % [mm]
    z = z_s*1e3; % [mm]
    nexttile(t1,iIm);
    imagesc(x,z,sono(:,:,1),1.5*[-1 1])
    colormap(sonomap)
    colorbar
    axis equal
    axis tight
    title(['f_v=',num2str(VibFreqArray(iIm)), 'Hz'])
    %%%%%%%%%%%%%% PLOT PATTERN %%%%%%%%%%%%%%

    %%%%%%%%%%%%%% PLOT B-MODE %%%%%%%%%%%%%%
    x = Properties.Width_B*1e3; % [mm]
    z = Properties.Depth_B*1e3; % [mm]
    nexttile(t2,iIm);
    imagesc(x,z,Properties.Bmode, caxis_bmode)
    colormap gray
%     colorbar
    axis equal
    axis tight
    title(['f_v=',num2str(VibFreqArray(iIm)), 'Hz'])
    %%%%%%%%%%%%%% PLOT B-MODE %%%%%%%%%%%%%%

%     save([sonoPath,'/',num2str(iIm),'.mat'],"sono","Properties");
end 

end

save_all_figures_to_directory(outiniPath)
%% --------------------------------------------------------------------- %%
%% 2-DIMENSIONAL TOTAL VARIATION
Nim = 9; % Number of channels
% swsPath = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\TV\sws';
ParamsTV.mu = 5;
ParamsTV.tol = 1e-1;
ParamsTV.isotropic = true;

RW_Param.dual = boolean(1); RW_Param.k = 1;
RW_Param.N_window = 20; RW_Param.beta = 1/100000;
RW_Param.tolerance = 1/1000;RW_Param.operator = 'G';
RW_Param.alpha=2;

% fig1 = figure('Position',[100 100 1000 500]);
% t1 = tiledlayout(fig1,3,3);
% fig2 = figure('Position',[100 100 1000 500]);
% t2 = tiledlayout(fig2,3,3);
% fig3 = figure('Position',[100 100 1000 500]);
% t3 = tiledlayout(fig3,3,3);

% Looping
% Nacq = input('Enter # of acquisitions: ');
Nacq = 3; % Number of acquisitions
% Nim = input('Enter # of freq. channels: ');
Nim = 20; % Number of channels
VibFreqArray = 120:20:500; % [Hz]

for iacq = 1:Nacq

fig1 = figure(10*iacq+2);
set(10*iacq+2, 'Position',[100 50 1500 800]);
t1 = tiledlayout(fig1,4,5);
sgtitle(['\bf SWS R-WAVE Method Acq #', num2str(iacq)], 'FontSize', font );

fig2 = figure(10*iacq+3);
set(10*iacq+3, 'Position',[100 50 1500 800]);
t2 = tiledlayout(fig2,4,5);
sgtitle(['\bf SWS TV Method Acq #', num2str(iacq)], 'FontSize', font );


for iIm = 1:Nim
%     load([sonoPath,'\',num2str(iIm),'.mat'])
    %load([swsPath,'\',num2str(iIm),'.mat']);


    load([dataDir, 'MatlabProcessed/Image', num2str(10*iacq+iIm),'/sono_filt.mat']);

    % Subsampling
    R = 2; % DECIMATION FACTOR
    [Nz,Nx,Nt] = size(sono_filt);
    sonoNew = zeros([ceil(Nz/R),Nx,Nt]);
    [b,a] = butter(4,1/R);
    for ix = 1:Nx
        for it = 1:Nt
            signal = filtfilt(b,a,sono_filt(:,ix,it));
            sonoNew(:,ix,it) = signal(1:R:end);
        end
    end
    
    % Selecting ROI
    zNew = z_s(1:2:end);
    izROI = zNew>z0 & zNew<zf;
    sono = sonoNew(izROI,:,:);
    Properties.Depth_S = zNew(izROI);
    izBmode = z_b >z0 & z_b <zf;
    Properties.Bmode = Bmode(izBmode,:);
    Properties.Depth_B = z_b(izBmode);
    Properties.Width_B = x_b;
    
    x = x_s*1e3; % [mm]
    z = z_s*1e3; % [mm]


    Properties.VibFreq = VibFreqArray(iIm);
    Properties.pitch = 3.0800e-04; % [m]

    tic
    swsRW = rWave2(sono,Properties,RW_Param);
    t = toc;
    fprintf('Exec. time for R-W: %f\n',t)

    fprintf("\nVibration Frequency = %d Hz\n",Properties.VibFreq);
    [swsTV,C,lin_system] = calculateSWSTV(sono,Properties,ParamsTV);
%     swsTV=1; C=1; linSystem=1;
    fprintf("Number of iterations: %d\n",length(C));
    

    
    tic
%     swsCWT = process_CWT(sono,Properties,[1 16]); % nnot yet
    swsCWT = 1;
    t = toc;
    fprintf('Exec. time for CWT: %f\n',t)
    
%     save([swsPath,'/',num2str(iIm),'.mat'],...
%         'swsTV','C','swsRW','swsCWT','Properties');

    % Selecting ROI
%     x = Properties.Width_S*1000;
%     z = Properties.Depth_S*1000;
    SWS_im_range = [2 6];
    
    % PLOT SWS R-WAVE
    nexttile(t1,iIm);
    imagesc(x,z,swsRW,SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['\alpha=',num2str(RW_Param.alpha,2),...
        ', f_v=',num2str(Properties.VibFreq), 'Hz'])
    ax = gca; ax.FontSize = 12;

    % PLOT SWS TV
    nexttile(t2,iIm);
    imagesc(x,z,swsTV,SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
%     title(['\mu=',num2str(ParamsTV.mu,2),...
%         ', f_v=',num2str(Properties.VibFreq),'Hz'])
    title([...
        'f_v=',num2str(Properties.VibFreq),'Hz'])
    ax = gca; ax.FontSize = 12;
 
    % PLOT SWS CWT
%     nexttile(t3,iIm);
%     swsCwtIm = squeeze(mean(swsCWT,3));
%     imagesc(x,z,swsCwtIm,SWS_im_range);
%     colormap turbo
%     colorbar
%     axis equal
%     xlim([x(1) x(end)]), xlabel('x [mm]')
%     ylim([z(1) z(end)]), ylabel('z [mm]')
%     title(['SWS from CWT, f_v=',num2str(Properties.VibFreq)])
%     ax = gca; ax.FontSize = 12;
end

end
save_all_figures_to_directory(outiniPath);

%% Selecting ROI
x0inc = 15; z0 = 11; L = 11; x0back = 1.5;
figure, 
%imagesc(x,z,swsTV,SWS_im_range); colormap turbo
imagesc(Properties.Width_B*1e3,Properties.Depth_B*1e3,...
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
    load([swsPath,'/',num2str(iIm),'.mat']);

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
% selectedImages = [3,4,5,7,8,9];
selectedImages = 1:9;
Nim = length(selectedImages);
for iIm = 1:Nim
    % Loading
    load([sonoPath,'\',num2str(selectedImages(iIm)),'.mat'])
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
        A = [A,                 sparse(Htot,Wim);
            sparse(Him,Wtot),   Aim];
        B = [B;Bim];
    end
    t = toc;
    fprintf('Exec. time for merging the systems: %f\n',t)
end

save([swsPath,'\TNV_all.mat'],'A','B');
%% Algorithm
% larger tau increases num of iterations and makes the results less
% dependendent on lambda
% smaller tau increases agreement between channels
lambda = 5; % +-50%
tau = 0.5;   % +-50%
maxIter = 500;
tol = 3e-4;
stableIter = 50;
Nim = length(selectedImages);
tic
[u, cost, error, fide, regul] = pdo_inv_tnv(B, Nz, Nx, A, ...
    lambda, tau, maxIter, tol, Nim, stableIter);
swsTNV = reshape(u,Nz,Nx,Nim);
toc
save([swsPath,'\TNV2.mat'],'swsTNV','Properties');


%% Results
load([swsPath,'\TNV2.mat'],'swsTNV');
x = Properties.Width_S * 1000;
z = Properties.Depth_S * 1000;
SWS_im_range = [2,6];
VibFreqArray = 200:20:360;

% figure('Position',[100 100 1200 600]);
figure,
t = tiledlayout(3,3);
% t = tiledlayout(2,3);
for iIm = 1:Nim
    nexttile
    imagesc(x,z,squeeze(swsTNV(:,:,iIm)),SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['TNV Inv F_v = ',num2str(VibFreqArray(selectedImages(iIm))),' Hz'])
    ax = gca; ax.FontSize = 12;
end
sgtitle(t,['\bfSWS from TNV, \lambda=',num2str(lambda,2),...
    ', \tau=',num2str(tau,2)])

%% DENOISING FEBRURARY
% figure, 
% montage(swsTNV), colormap jet


%%%%%%%%%%%%%%%%%     OPTIMIZED FEBRURARY 2023 %%%%%%%%%%%%%%%%%
    reshapedSWS = reshape(swsTNV, [], size(swsTNV, 3));
    meanValues = mean(reshapedSWS, 1, 'omitnan');
    stdValues = std(reshapedSWS, 0, 1, 'omitnan');
    
    ratios = abs(meanValues) ./ stdValues;
    SWS.SNR_ratios = ratios;

    figure, plot(VibFreqArray, SWS.SNR_ratios , 'o-'), grid on; 
    title('SNR'), xlabel('Freq [Hz]')
%%%%%%%%%%%%%%%%%     OPTIMIZED FEBRURARY 2023 %%%%%%%%%%%%%%%%%

lambda = 2; % +-50%
tau = 0.5;   % +-50%
maxIter = 500;
tol = 3e-4;
stableIter = 50;
snr_weight = 0;
if (snr_weight) weightEstimators =  SWS.SNR_ratios;
else weightEstimators =  ones(size(SWS.SNR_ratios)); end

[x, cost, error, fid, reg] = pdo_den_wtnv(swsTNV, lambda, tau, maxIter, tol, stableIter, weightEstimators);
sws_denTNV = x;

% Results DENOISING
load([swsPath,'\TNV2.mat'],'swsTNV');
x = Properties.Width_S * 1000;
z = Properties.Depth_S * 1000;
SWS_im_range = [2,6];
VibFreqArray = 200:20:360;

% figure('Position',[100 100 1200 600]);
figure,
t = tiledlayout(3,3);
% t = tiledlayout(2,3);
for iIm = 1:Nim
    nexttile
    imagesc(x,z,squeeze(sws_denTNV(:,:,iIm)),SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['TNV Den Ch_{SNR}', num2str(snr_weight),' F_v = ',num2str(VibFreqArray(selectedImages(iIm))),' Hz'])
    ax = gca; ax.FontSize = 12;
end
sgtitle(t,['\bfSWS from TNVDen, \lambda=',num2str(lambda,2),...
    ', \tau=',num2str(tau,2)])


%% Results
load([swsPath,'\TNV2.mat'],'swsTNV');
x = Properties.Width_S * 1000;
z = Properties.Depth_S * 1000;
SWS_im_range = [2,6];
VibFreqArray = 200:20:360;

% figure('Position',[100 100 1200 600]);
figure,
t = tiledlayout(3,3);
% t = tiledlayout(2,3);
for iIm = 1:Nim
    nexttile
    imagesc(x,z,squeeze(swsTNV(:,:,iIm)),SWS_im_range);
    colormap turbo
    colorbar
    axis equal
    xlim([x(1) x(end)]), xlabel('x [mm]')
    ylim([z(1) z(end)]), ylabel('z [mm]')
    title(['F_v = ',num2str(VibFreqArray(selectedImages(iIm))),' Hz'])
    ax = gca; ax.FontSize = 12;
end
sgtitle(t,['\bfSWS from TNV, \lambda=',num2str(lambda,2),...
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

%% DENOISING 

%% Selecting ROI
% swsPath = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\TV\sws';
load([swsPath,'\1.mat']);
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
[X,Z] = meshgrid(Properties.Width_S*1e3,Properties.Depth_S*1e3);
maskInc = (X>x0inc & X<x0inc+L & Z>z0 & Z<z0+L);
maskBack = (X>x0back & X<x0back+L & Z>z0 & Z<z0+L);

cnrTV = NaN(1,Nim);
cnrRW = NaN(1,Nim);
cnrTNV = NaN(1,Nim);

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
end
%% Comparing CNR
VibFreqArray = 200:20:360;
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
    x = Properties.Width_S*1000;
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
    zi = Properties.Depth_S(iz)*1000;
    title(['Lateral profile at z=',num2str(zi,2),'mm, f_v=',...
        num2str(freqArr(iIm)),'Hz'])
end
%% Comparing axial profile
ix = 80;
freqArr = 200:20:360;
for iIm = 1:Nim
    load([swsPath,'\',num2str(iIm),'.mat']);
    z = Properties.Depth_S*1000;
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
    xi = Properties.Width_S(ix)*1000;
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
    x = Properties.Width_S(ix0:ixf-1)*1000;
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
    zi = Properties.Depth_S(iz)*1000;
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
ylim([2 14])

%% -------------------------------------------------------------------- %%
%% USEFUL FUNCTIONS
function [swsTV,C, linSystem] = calculateSWSTV(sono,Properties,ParamsTV)
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
linSystem.A = A;
linSystem.b = B;
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