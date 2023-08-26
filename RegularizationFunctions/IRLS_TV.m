% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function [u,C] = IRLS_TV(b,A,mu,M,N,tol,mask,isotropic,colMajor)
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
%       G               vector containing cost function for each iteration
%
% Author: Andres Leonel Coila
% Modified by Sebastian Merino

AtA = A'*A; % This takes A LOT OF TIME
Atb = A'*b;
    
[u,~] = cgs2(AtA,Atb,1e-6,20);
C(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask,isotropic,colMajor);

[Dx,Dy] = diffOperator(M,N,colMajor);
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
        
        C(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask,isotropic,colMajor);
        error = abs(C(ite_irls+1) - C(ite_irls))/C(ite_irls+1);
    end

else
    while error > tol && ite_irls < 1000        
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
        
        C(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask,isotropic,colMajor);
        error = abs(C(ite_irls+1) - C(ite_irls))/C(ite_irls+1);
    end
end

end

function [Dx,Dy] = diffOperator(M,N,colMajor)
% Returns differential operators
% Inputs: 
%       M,N             image size
%       colMajor        true if the matrix is stored col-major, false if
%                       stored row-major
%  
% Outputs:
%       Dx,Dy           Difference operators
% Author: Sebastian Merino
if ~colMajor
    foo = M;
    M = N;
    N = foo;
end

G = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
G(:,end) = [];
G(M,M) = 0;
Dx = kron(speye(N),G);

G = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
G(:,end) = [];
G(N,N) = 0;
Dy = kron(G,speye(M));

end