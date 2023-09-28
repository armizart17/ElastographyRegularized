% Script para procesar data de adquisision WIDE BEAM 27-09
clear, clc
addpath('./RegularizationFunctions/',"./ElastographyFunctions")

baseDir = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\27_09_23';
rawDir = [baseDir,'\raw'];
rawFiles = dir([rawDir,'\*.mat']);
Nim = length(rawFiles);

iqDir = [baseDir,'\iq'];
pvDir = [baseDir,'\pv'];
swsDir = [baseDir,'\sws'];

% Colormap
load('MyColormaps.mat')
%% Processing raw files
if ~exist(iqDir,'dir'), mkdir(iqDir); end
for iAcq = 1:Nim
    file = rawFiles(iAcq);
    load([file.folder,'\',file.name]);
    
    % Acquisition parameters
    dinf.vibf = str2double(file.name(7:9));
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
    [Nz,Nx,Nens] = size(IQ);
    
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
    maskIQ = logical(repmat(mask,[1,1,Nens]));
    IQ = reshape(IQ(maskIQ),length(dinf.z),length(dinf.x),Nens);
    IQ = IQ(1:end-1,:,:,:);
    dinf.z = dinf.z(1:end-1);
    
    %Plotting B-mode images
    figure, tiledlayout(1,2);
    nexttile,
    x = dinf.xBmode*1000;
    z = dinf.zBmode*1000;
    imagesc(x,z,dinf.Bmode)
    title('Bmode')
    axis equal
    axis tight
    colormap gray

    nexttile
    x = dinf.x*1000;
    z = dinf.z*1000;
    imagesc(x,z,db(IQ(:,:,1,1)))
    title('Doppler area')
    axis equal
    axis tight
    colormap gray
    save([iqDir,'\',num2str(iAcq),'.mat'],'IQ','dinf')
end

%% Particle velocity estimator
if ~exist(pvDir,'dir'), mkdir(pvDir); end

for iAcq = 1:Nim
    load([iqDir,'\',num2str(iAcq),'.mat'])

    tic
    [pv,dinf] = pv_cal(IQ,dinf,1);
    toc
    
    % Displaying pv
    [Nz,Nx,Nt] = size(pv);
    amp = 3*std(pv(:));
    figure,
    for it = 1:10
        imagesc(pv(:,:,it), amp*[-1,1])
        pause(1/5)
    end
    
    % Point in time
    ix = 40;
    iz = 50;
    t = (0:Nt-1)/dinf.PRF;
    pvPoint = squeeze(pvFilt(iz,ix,:));
    Nf = 128;
    pvFT = fft(pvPoint,Nf);
    f = (0:Nf-1)/Nf*dinf.PRF;
    
    figure,tiledlayout(2,1)
    nexttile
    plot(t,pvPoint)
    xlabel('Time [s]')
    ylabel('PV [m/s]')
    title('Time domain')
    
    nexttile
    plot(f,abs(pvFT))
    xlabel('Frequency [Hz]')
    title('Freq. domain')

    save([pvDir,'\',num2str(iAcq),'.mat'],'pv','dinf')
end


%% Filtering
if ~exist(swsDir,'dir'), mkdir(swsDir); end

for iAcq = 1:Nim
    load([pvDir,'\',num2str(iAcq),'.mat'])

    sd = [2,7];
    [pvFilt,filtMask] = directionalFilter(pv,dinf.PRF,dinf.dx,...
        dinf.vibf,dinf.vibf,sd,'left');
    amp = 3*std(pvFilt(:));
    figure,
    for it = 1:10
        imagesc(dinf.x, dinf.z, pvFilt(:,:,it), amp*[-1,1])
        pause(1/10)
    end
    
    % CWT
    swsRange = [3 7];
    fLim = dinf.vibf./[swsRange(2) swsRange(1)]; % SWS range
    fMax = freq_CWT(pv,1/dinf.dx,fLim);
    
    swsCWT = dinf.vibf ./ fMax;
    
    x = dinf.x*1000;
    z = dinf.z*1000;
    swsImRange = [3 7];
    figure,
    imagesc(x,z,swsCWT(:,:,1),swsImRange);
    colormap turbo
    colorbar
    axis equal
    axis tight
    xlabel('x [mm]'), ylabel('z [mm]')
    ax = gca; ax.FontSize = 12;

    save([swsDir,'\',num2str(iAcq),'.mat'],'swsCWT','dinf')

end
