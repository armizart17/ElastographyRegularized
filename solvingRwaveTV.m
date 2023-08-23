clear, clc
addpath('./RegularizationFunctions/',"./ElastographyFunctions")
load data360Hz.mat
%% Testing function on site
b = B;
mu = 1;
[M,N,~] = size(sonoSub);
tol = 1e-1;
mask = ones(M*N,1);
isotropic = true;
colMajor = false;

AtA = A'*A; % This takes A LOT OF TIME
Atb = A'*b;
    
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
        error = abs(G(ite_irls+1) - G(ite_irls));
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

