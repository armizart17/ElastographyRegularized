% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_TV(b,A,mu,M,N,tol,mask,isotropic,colMajor)
% Returns the Total Variation
% Inputs: 
%       b               vector containing measurements
%       A               matrix for linear system of eq
%       mu              regularization parameter
%       M,N             image size of u
%       tol             tolerance for relative error
%       mask            binary mask of analyzed region
%       isotropic       true if the TV is isotropic, false if anisotropic
%       colMajor        true if the matrix is stored col-major, false if
%                       stored row-major
%  
% Outputs:
%       u               vector of image samples, size MN
%
% Author: Andres Leonel Coila
% Modified by Sebastian Merino

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
eps = 0.01; % 0.01 for isotropic is ok
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
        %fprintf("IRLS iteration %d\n",ite_irls);
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