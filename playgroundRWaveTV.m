%% Loading filtered sono signal
clear, clc
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
figure,
for it = 1:Nt
    imagesc(x,z,sonoFilt(:,:,it),1.5*[-1 1])
    colorbar
    axis equal
    axis tight
    pause(1/Properties.FrameRate)
end
%% Subsampling
sonoNew = zeros([Nz/2,Nx,Nt]);
for ix = 1:Nx
    for it = 1:Nt
        sonoNew(:,ix,it) = decimate(sonoFilt(:,ix,it),2);
    end
end

%% Playing subsampled video
figure,
for it = 1:Nt
    imagesc(x,z,sonoNew(:,:,it),1.5*[-1 1])
    colorbar
    axis equal
    axis tight
    pause(1/Properties.FrameRate)
end

%%

sonoNew2 = zeros([Nz/2,Nx,floor(Nt/2)]);
for ix = 1:Nx
    for it = 1:Nt
        sonoNew2(:,ix,it) = decimate(sonoFilt(:,ix,it),2);
    end
end


