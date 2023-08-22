M = 10;
N = 6;
D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;

figure, imagesc(D)
colorbar
axis equal
axis tight
%% Differential operator in x
Dx = kron(speye(N),D);
figure, imagesc(Dx)
colorbar
axis equal
axis tight

%% Differential operator in y
D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));
figure, imagesc(Dy)
colorbar
axis equal
axis tight
%% D matrix
D = [Dx' Dy']';
figure, imagesc(D)
colorbar
axis equal
axis tight

%% Sample matrix X
X = ones(M,N);
X(1:floor(M/2),1:floor(N/2)) = 0;
figure, imagesc(X)
colorbar
axis equal
axis tight
%% Difference matrices and differential operators
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

figure, 
subplot(221), 
imagesc(reshape(Dh,[M,N]))
colorbar
axis equal
axis tight
title('Horizontal diff')

subplot(222), 
imagesc(reshape(Dv,[M,N]))
colorbar
axis equal
axis tight
title('Vertical diff')

subplot(223), 
imagesc(reshape(Dx*X(:),[M,N]))
colorbar
axis equal
axis tight
title('Dx operator')

subplot(224), 
imagesc(reshape(Dy*X(:),[M,N]))
colorbar
axis equal
axis tight
title('Dy operator')

%% Inverse of gradients
P = Dh.^2 + Dv.^2;
eps = 0.01;
P = 2*sqrt(P.^2 + eps^2);
P = P.^(-0.5);
figure, imagesc(P)
colorbar
axis equal
axis tight
%% W and omega matrices
P = P(:);
omega = speye(M*N);
omega = spdiags(P,0,omega);
W = kron(speye(2),omega);
figure, imagesc(W)
colorbar
axis equal
axis tight

%% Solving
A = randi(10,[1000,length(X(:))]);
b = randi(10,[1000,1]);
mu = 1;
size(A'*A + mu*D'*W*D)
size(X(:))
size(A'*b)