% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_TV(b,A,mu,M,N,tol,mask,isotropic)

[u,~] = cgs2(A'*A,A'*b,1e-6,20);

G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask,isotropic);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dx = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

if (isotropic)
    while error > tol
        
        X = reshape(u,M,N);
        ite_irls = ite_irls + 1;
        Dh = [diff(X,[],1);zeros(1,N)];
        Dv = [diff(X,[],2) zeros(M,1)];
        
        P = Dh.^2 + Dv.^2;
        eps = 0.01;
        P = 2*sqrt(P.^2 + eps^2);
        P = P.^(-0.5);
        P = P(:).*mask;
        omega = speye(M*N);
        omega = spdiags(P,0,omega);
        W = kron(speye(2),omega);
        
        AtA = A'*A;
        [u] = cgs2( AtA + mu*D'*W*D, A'*b, 1e-6 , 20, u );
        
        G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask,isotropic);
        error = abs(G(ite_irls+1) - G(ite_irls));
        
    end

else
    % Iteratively reweighted least squares (anisotropic)
    while error > tol && ite_irls < 20
        
        X = reshape(u,M,N);
        ite_irls = ite_irls + 1;
        Dh = [diff(X,[],1);zeros(1,N)];
        Dv = [diff(X,[],2) zeros(M,1)];
        
        eps = 0.1;
        
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
        
        AtA = A'*A;
        [u] = cgs2( AtA + mu*Dx'*Wx*Dx + mu*Dy'*Wy*Dy , A'*b, 1e-6 , 20, u );
        
        G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask,isotropic);
        error = abs(G(ite_irls+1) - G(ite_irls));
        
    end
end
%figure(909); plot(1:length(G),G);

end