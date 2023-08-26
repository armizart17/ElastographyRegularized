clear, clc
addpath('./RegularizationFunctions/',"./ElastographyFunctions")
load data360Hz.mat
%% Testing function on site
b = B;
mu = 1;
[M,N,~] = size(sonoSub);
tol = 1e-3;
mask = ones(M*N,1);
isotropic = true;
colMajor = false;
tic
AtA = A'*A; % This takes A LOT OF TIME
Atb = A'*b;
toc
[u,~] = cgs2(AtA,Atb,1e-6,20);

G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask,isotropic,colMajor);

if colMajor
    D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
    D(:,end) = [];
    D(M,M) = 0;
    Dx = kron(speye(N),D);
    
    D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
    D(:,end) = [];
    D(N,N) = 0;
    Dy = kron(D,speye(M));
else
    D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
    D(:,end) = [];
    D(N,N) = 0;
    Dx = kron(speye(M),D);
    
    D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
    D(:,end) = [];
    D(M,M) = 0;
    Dy = kron(D,speye(N));
end
D = [Dx' Dy']';

ite_irls = 0;
error = 1;
eps = 0.1; % 0.01 for isotropic is ok

if (isotropic)
    while error > tol        
        ite_irls = ite_irls + 1;
        fprintf("IRLS iteration %d\n",ite_irls);

        Dh = Dx*u;
        Dv = Dy*u;
        P = sqrt(Dh.^2 + Dv.^2 + eps^2);
        P = P.^(-0.5);
        P = P(:).*mask;
        
        omega = speye(M*N);
        omega = spdiags(P,0,omega);
        
        W = kron(speye(2),omega);
        [u,~] = cgs2( AtA + mu*D'*W*D, Atb, tol , 20, u );
        
        G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask,isotropic,colMajor);
        error = abs(G(ite_irls+1) - G(ite_irls))/G(ite_irls+1);
    end

else
    while error > tol && ite_irls < 30        
        ite_irls = ite_irls + 1;
        fprintf("IRLS iteration %d\n",ite_irls);
        Dh = Dx*u;
        Dv = Dy*u;

        Px = abs(Dh + eps);
        Px = 1./Px;
        Px = Px(:).*mask;
        omega = speye(M*N);
        omega = spdiags(Px,0,omega);
        Wx = kron(speye(1),omega);
        
        Py = abs(Dv + eps);
        Py = 1./Py;
        Py = Py(:).*mask;
        omega = speye(M*N);
        omega = spdiags(Py,0,omega);
        Wy = kron(speye(1),omega);
        
        [u] = cgs2( AtA + mu*Dx'*Wx*Dx + mu*Dy'*Wy*Dy , Atb, tol, 20, u );
        
        G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask,isotropic,colMajor);
        error = abs(G(ite_irls+1) - G(ite_irls));
    end
end
%% Error per iteration
figure('Position', [200 200 500 500]);
plot(G)
xlabel("# of iterations")
ylabel("Error")
%%  Displaying and comparing results
baseDir = 'C:\Users\sebas\Documents\MATLAB\Elastography';
sonoPath = [baseDir,'\heterogeneo_sono'];

load([sonoPath, '\Image9\sono.mat'],'Properties');
dx = 3.0800e-04;
dz = (Properties.Depth_S(2) - Properties.Depth_S(1))*2;
Properties.pitch = dx;

swsTV = reshape(u,N,M)';
SWS_im_range = [2 5.5];
x = 1000*(1:N)*dx;
z = 1000*(20:140)*dz;
figure('Position', [200 200 500 500]);
subplot(211),
imagesc(x,z,swsTV,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title('SWS from TV')
ax = gca; ax.FontSize = 12;

RW_Param.dual = boolean(1); RW_Param.k = 1;
RW_Param.N_window = 20; RW_Param.beta = 1/100000;
RW_Param.tolerance = 1/1000;RW_Param.operator = 'G';
RW_Param.alpha=2;
tic
swsRWave = rWave2(sonoSub,Properties,RW_Param);
toc
subplot(212),
imagesc(x,z,swsRWave,SWS_im_range);
colormap turbo
colorbar
axis equal
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title('SWS from Tikhonov Reg')
ax = gca; ax.FontSize = 12;

%% Looping for reg parameter
dataPath = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\TV\muData\';
b = B;
muArray = logspace(-2,2,20);
[M,N,~] = size(sonoSub);
tol = 1e-1;
mask = ones(M*N,1);
isotropic = true;
colMajor = false;

for iMu = 1:length(muArray)
    mu = muArray(iMu);
    [u,G] = IRLS_TV(b,A,mu,M,N,tol,mask,isotropic,colMajor);
    swsTV = reshape(u,N,M)';
    save([dataPath,num2str(iMu),'.mat'],'swsTV','G','mu');
end

%% Getting figures for each reg parameter
muArray = logspace(-2,2,20);
figurePath = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\TV\muFigures\';
for iMu = 1:20
    load([dataPath,num2str(iMu),'.mat']);
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
    saveas(gcf,[figurePath,num2str(iMu),'.png'])
end

%% Selecting ROI
load([dataPath,num2str(13),'.mat']);
x0inc = 16; z0 = 12; L = 11; x0back = 1.5;

imagesc(x,z,swsTV,SWS_im_range);
colormap turbo
colorbar
axis equal
hold on
rectangle('Position',[x0inc z0 L L], 'LineWidth',2),
rectangle('Position',[x0back z0 L L], 'LineWidth',2, 'EdgeColor','w'),
hold off
xlim([x(1) x(end)]), xlabel('x [mm]')
ylim([z(1) z(end)]), ylabel('z [mm]')
title(['SWS from TV, \mu=',num2str(mu,2)])
ax = gca; ax.FontSize = 12;

%% Calculating CNR
[X,Z] = meshgrid(x,z);
maskInc = (X>x0inc & X<x0inc+L & Z>z0 & Z<z0+L);
maskBack = (X>x0back & X<x0back+L & Z>z0 & Z<z0+L);
cnrArray = zeros(1,20);
meanInc = zeros(1,20);
meanBack = zeros(1,20);
stdInc = zeros(1,20);
stdBack = zeros(1,20);
for iMu = 1:20
    load([dataPath,num2str(iMu),'.mat']);
    swsInc = swsTV(maskInc);
    swsBack = swsTV(maskBack);
    cnrArray(iMu) = 2*(mean(swsBack) - mean(swsInc))^2 / ...
        (var(swsInc) + var(swsBack));
    meanInc(iMu) = mean(swsInc);
    meanBack(iMu) = mean(swsBack);
    stdInc(iMu) = std(swsInc);
    stdBack(iMu) = std(swsBack);  
end
figure,
semilogx(muArray,cnrArray,'o-')
ylabel('CNR'), xlabel('mu')
%% SWS values for each mu
figure,
errorbar(muArray,meanBack,stdBack, 'LineWidth',2)
hold on
errorbar(muArray,meanInc,stdInc, 'LineWidth',2)
hold off
ax = gca;
ax.XScale = "log";
legend({'Back','Inc'})
ylabel('SWS [m/s]'), xlabel('mu')

%% Testing anisotropic
mu = 10;
tol = 1e-3;
isotropic = false;
[u,G] = IRLS_TV(b,A,mu,M,N,tol,mask,isotropic,colMajor);
swsTV = reshape(u,N,M)';
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
saveas(gcf,[figurePath,num2str(iMu),'.png'])
