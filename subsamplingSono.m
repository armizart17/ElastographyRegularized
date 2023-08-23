%% Loading filtered sono signal
clear, clc
addpath('./RegularizationFunctions/',"./ElastographyFunctions")
baseDir = 'C:\Users\sebas\Documents\MATLAB\Elastography';
%baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\Datasets';
sonoPath = [baseDir,'\heterogeneo_sono'];
swsRange = [2,8];  
move = 'left';

% Pre-processing
load([sonoPath, '\Image9\sono.mat']);
Properties.pitch = 3.0800e-04;
[sonoFilt,~,~] = process_sono_data(sono,Properties,move,swsRange);

%% Playing video
[Nz,Nx,Nt] = size(sonoFilt);
x = Properties.Width_S * 1000;
z = Properties.Depth_S * 1000;
figure('Position', [200 200 300 300]);
for it = 1:5
    imagesc(x,z,sonoFilt(:,:,it),1.5*[-1 1])
    colorbar
    colormap(sonomap)
    title('Sonoelastography')
    xlabel('Lateral position [mm]')
    ylabel('Depth [mm]')
    axis equal
    axis tight
    pause(1/Properties.FrameRate)
end
%% Subsampling in depth
sonoNew = zeros([Nz/2,Nx,Nt]);
for ix = 1:Nx
    for it = 1:Nt
        sonoNew(:,ix,it) = decimate(sonoFilt(:,ix,it),2);
    end
end

%% Playing subsampled video
figure('Position', [200 200 300 300]);
for it = 1:5
    imagesc(x,z,sonoNew(:,:,it),1.5*[-1 1])
    colorbar
    colormap(sonomap)
    title('Sonoelastography subsampled')
    xlabel('Lateral position [mm]')
    ylabel('Depth [mm]')
    axis equal
    axis tight
    pause(1/Properties.FrameRate)
end

%% Subsampling in time
sonoNew2 = zeros([Nz/2,Nx,ceil(Nt/2)]);
for ix = 1:Nx
    for iz = 1:Nz/2
        sonoNew2(iz,ix,:) = decimate(squeeze(sonoNew(iz,ix,:)),2);
    end
end
%%
figure,
for it = 1:size(sonoNew2,3)
    imagesc(x,z,sonoNew2(:,:,it),1.5*[-1 1])
    colorbar
    axis equal
    axis tight
    pause(1/Properties.FrameRate)
end
%% Computing SWS images
RW_Param.dual = boolean(1); RW_Param.k = 1;
RW_Param.N_window = 20; RW_Param.beta = 1/100000;
RW_Param.tolerance = 1/1000;RW_Param.operator = 'G';
RW_Param.alpha=2;
tic
SWSsub = rWave2(sonoNew2,Properties,RW_Param);
toc
tic
SWS = rWave2(sonoFilt,Properties,RW_Param);
toc
%% Comparing results
SWS_im_range = [2 5.5];
x = 1000*Properties.Width_S;
z = 1000*Properties.Depth_S;
figure('Position', [200 200 800 350]);
subplot(121),
imagesc(x,z,SWS,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title('SWS from raw video')
ax = gca; ax.FontSize = 12;

subplot(122),
imagesc(x,z,SWSsub,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title('SWS from subsampled video')
ax = gca; ax.FontSize = 12;

%% Selecting ROI and cleaning workspace
sonoSub = sonoNew2(20:140,:,:);
clear sono sonoNew2 sonoNew sonoFilt

%% Generating A matrix
[Nz,Nx,Nt] = size(sonoSub);
%iz = floor(Nz/2);
B = [];
A = [];
for iz = 1:Nz
    Az = zeros(1,Nx); % Weight matrix for line iz
    n=1;
    for frame = 1:size(sonoSub,3)
        % Finding indices for each peak in each line (1x128)
        sonoLine = squeeze(sonoSub(iz,:,frame));
        [iPeaks] = peakfinder(sonoLine/max(sonoLine(:)),0);
        
        % Ignoring if peaks are in the start, the end, or there are none
        if size(iPeaks,2)<2, continue; end
        if(iPeaks(1)==1), iPeaks=iPeaks(2:end); end          
        if(iPeaks(end)==size(sonoSub,2)), iPeaks=iPeaks(1:end-1); end
        if size(iPeaks,2)<2, continue; end
    
        % SWS for each wavelength
        lamdaSamples = diff(iPeaks);
        swsLamda = lamdaSamples*2* Properties.pitch * Properties.VibFreq;
    
        % Extrapolation at the start (nearest neighbour)
        Az(n,1:iPeaks(1)) = 1/(iPeaks(1));
        B = [B; swsLamda(1)];
        n = n+1;
    
        % Estimation of weighted coefficients
        for i=1:size(swsLamda,2) % For each wavelength
            Az(n,iPeaks(i)+1:iPeaks(i+1))= 1/lamdaSamples(i);
            B = [B; swsLamda(i)]; 
            n = n+1;
        end
        
        % Extrapolation at the end
        Az(n,(iPeaks(end)+1):Nx) = 1/(Nx-(iPeaks(end)+1)+1); 
        B = [B; swsLamda(end)]; 
        n = n+1;
    end
    
    if RW_Param.dual
        % Find valleys too
        for frame=1:size(sonoSub,3)
            sonoLine=squeeze(-1*sonoSub(iz,:,frame));
    
            % Repeating the same procedure
            [iPeaks] = peakfinder(sonoLine/max(sonoLine(:)),0);
            if size(iPeaks,2)<2, continue; end
            if(iPeaks(1)==1), iPeaks=iPeaks(2:end); end
            if(iPeaks(end)==size(sonoSub,2)), iPeaks=iPeaks(1:end-1); end
            if size(iPeaks,2)<2, continue; end
            lamdaSamples = diff(iPeaks);
            swsLamda = lamdaSamples*2* Properties.pitch * Properties.VibFreq;
            Az(n,1:iPeaks(1)) = 1/(iPeaks(1));
            B = [B; swsLamda(1)];
            n = n+1;
            for i=1:size(swsLamda,2) % For each wavelength
                Az(n,iPeaks(i)+1:iPeaks(i+1))= 1/lamdaSamples(i);
                B = [B; swsLamda(i)]; 
                n = n+1;
            end
            Az(n,(iPeaks(end)+1):Nx) = 1/(Nx-(iPeaks(end)+1)+1);
            B = [B; swsLamda(end)];
            n = n+1;
        end
    end
    %Nw = size(Az,1);    % Number of wavelengths
    if iz == 1
        A = Az;
    else
        A = [A zeros(size(A,1),Nx);zeros(size(Az,1),Nx*(iz-1)),Az];
    end
end
%% Plotting matrices
figure('Position', [100 100 500 700]), 
imagesc(A)
set(gca,'cLim',[0.04 0.08])
title('Weight matrix A')
xlabel('Lateral sample')
ylabel('# of wavelength')
colorbar
xlim([1,256])
ylim([1,500])

%% plotting sono Sub
[Nz,Nx,Nt] = size(sonoSub);
dx = 3.0800e-04;
dz = (Properties.Depth_S(2) - Properties.Depth_S(1))*2;
x = 1000*(1:Nx)*dx;
z = 1000*(20:140)*dz;

figure('Position', [100 100 600 400]), 
imagesc(x,z,sonoSub(:,:,1))
%set(gca,'cLim',[0.04 0.08])
title('Sonoelastography ROI')
xlabel('Lateral position [mm]')
ylabel('Depth [mm]')
colorbar
colormap(sonomap)
axis equal
axis tight
%xlim([1,256])
%ylim([1,500])


%%
save data360Hz.mat A B sonoSub -v7.3
