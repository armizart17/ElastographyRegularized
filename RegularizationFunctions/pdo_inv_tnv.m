function [x, cost, error, fide, regul] = pdo_inv_tnv(Y, A, lambda, tau, maxIter, tol, numberEstimators, stableIter)
% function [x, cost, error, fide, regul] = pdo_inv_tnv(Y, A, lambda, tau, maxIter, tol, numberEstimators, stableIter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Primal DualSpliting Optimization Method for linear inverse problem A*x = y 
% Cost function to minimize: arg min || Ax - y ||_2^2 + TNV (x)
% This was programmed to solve linear problems for images, i.e 2D or even 2D multichannel images
% INPUTS
%	- Y		: Observed data (Array form)
%	- A		: Linear operator (Array form, matrix)
% 	- lambda: regularization parameter (try powers of 10 for linear inverse problem)
% 	- tau	: second hyperparameter for proximal operators, apply grid search i.e. 0.1 0.01, 0.001
%	- maxIter: maximum number of iterations, set 10000
% 	- tol	: tolerance error for cost function
%   -numberEstimators: Defines the # of channels of the ouput image (i.e 2 for TNV-SLD)
% 	-stableIter: If you plot cost functions, depending on the problem the cost function has atypical values for 50-100 iterations
%				this could cause the algorithm to stop regularization before time, set this 50-100

% OUTPUTS
% 	- x: output data, it returns as an array, could be a 2D multichannel image
% 	- cost: cost function vector according to # iterations
%   - error: We are minimizing the cost function so we analyze the first derivative to set where the cost function is minimized
%      when error are small, cost function has achieved is local minimum.
%   - fide: fidelization term vector according to # iterations
%   - regul : regularization vector according to # iterations


% AUTHOR: Edmundo Arom Miranda
% VERSION: pdo_inv_tnv v1.0 for Sebastian Merino 

% REFERENCE:
% Based on 
% Condat, L. (2013). A primal–dual splitting method for convex optimization involving 
% Lipschitzian, proximable and linear composite terms. Journal of optimization theory and applications, 158(2), 460-479.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	rho = 1.99;		% relaxation parameter, in [1,2]
	sigma = 1/tau/8; % proximal parameter

    global H W C numEstim
    [H,W,C]=size(Y);
    numEstim = numberEstimators;

    y = Y(:);
    At = A';
    AtA = At*A;

    global inv_izq
    [m n] = size(A);
    
    izq = sparse(eye(n)) + tau*(AtA);
    inv_izq = ( izq )^-1;
    
	% INITIALIZATION OPTIONS TO HELP THE REGULARIZATION
    x_ini = (A'*A)\(A'*y);% QR solver 
%     x_ini = cgs(A'*A, A'*y);    % CGS solver 
%     x_ini = zeros([H, W, numEstim]); % zeros initialization
    
    
%   FAST IMPLEMENTATION FOR CC channels
    x2 = reshape(x_ini, [H,W,numEstim]);
   
    x2(:,:,1) = Y( :,:,floor(C/2) ); % 1st estimator

    u2 = zeros([size(x2) 2]);

	ee = 1;	
    error(1) = 1;
    iter = 1;
    cost(1) = 100;
    fide(1) = 1;
    regul(1) = 1;
   
    disp('Itn   |    Cost   |   Delta(Cost)');
    disp('----------------------------------');
    

    while (iter < maxIter) && (ee > tol)    

        x = prox_tau_f_generic(x2-tau*opDadj(u2),Y,tau, A, AtA, At);
		u = prox_sigma_g_conj(u2+sigma*opD(2*x-x2), lambda);
        
		x2 = x2 + rho*(x-x2);
		u2 = u2 + rho*(u-u2);
       
        
        fide(iter+1) = 0.5*sum(sum(sum((reshape(A*x(:), size(Y)) - Y).^2)));
        regul(iter+1) = +lambda*TVnuclear_EMZ(x);
        
        cost(iter+1) = fide(iter+1) + regul(iter+1);
%         cost(iter+1) = 0.5*sum(sum(sum((reshape(A*x(:), size(Y)) - Y).^2)))+lambda*TVnuclear_EMZ(x);
        
        ee = abs(cost(iter+1) - cost(iter));     
		
        % NO NORMALIZADED
%         error(iter+1) = ee;
        
        % NORMALIZED
        ee = ee / cost(iter+1);
		
        error(iter+1) = ee;
        
        if (iter < stableIter) % en las primeras 100 iteraciones inestable 
            ee = 1; % poner error grande (cualquier valor) para que no pare iteracion            
        end
    
        
		if mod(iter,5)==0 % to print every 5 iterations
%             primalcost = 0.5* norm(A*x(:) - y,2).^2 + lambda*TVnuclear_EMZ(x);         
            fprintf('%4d  | %f | %e\n',iter, cost(iter+1), error(iter+1));

        end
        iter = iter+1;
	end
return

function [u] = prox_tau_f(v, B, tau)
    
    u =  (v+tau*B)/(1+tau);

return

function [u_3d] = prox_tau_f_generic(v, B, tau, A, AtA, At)

    global H W C numEstim;
    global inv_izq
    [m n] = size(A);
    
    der = v(:)+tau*At*B(:);
    
    u = ( inv_izq )* ( der );
    
    u_3d = reshape(u, [H W numEstim]); % reshape as image
    
return
    
function [u] = prox_sigma_g_conj(v, lambda)

    u = v - prox_g(v,lambda);
return

function [u] = opD(v)
    % Description: Dy(:,:,1) Dx (:,:,:,2), 3rd dimension is channel
    [H,W,C]=size(v);
    % usually u = [Dx(v), Dy(v)]; for each channel use for
    
%     jacobian_vImage = DxDy(vectorial_image(:,:,1));
% 
%     for i = 2:size(vectorial_image,3)
%         v = vectorial_image(:,:,i);
%         jacobian_vImage = [jacobian_vImage; DxDy(v)];
%     end
    
    % one line
%    u = cat(4,[diff(v,1,1);zeros(1,W,3)],[diff(v,1,2) zeros(H,1,3)]);
    u = cat(4,[diff(v,1,1);zeros(1,W,C)],[diff(v,1,2) zeros(H,1,C)]);
return

function [u] = opDadj(v)
    % u = Dy Dx
    u = -[v(1,:,:,1);diff(v(:,:,:,1),1,1)]-[v(:,1,:,2) diff(v(:,:,:,2),1,2)];	
return
    
function [u] = opDadj_adap(v, factor) %% ADAPTATIVO
    % u = Dy Dx
%     u = -[v(1,:,:,1);diff(v(:,:,:,1),1,1)]-[v(:,1,:,2) diff(v(:,:,:,2),1,2)];
    
    u = factor .* (-[v(1,:,:,1);diff(v(:,:,:,1),1,1)]- ...
        [v(:,1,:,2) diff(v(:,:,:,2),1,2)]);
return

function [u] = TVnuclear_EMZ(v)

    % u = norm(eig(opD(v)),1); % oversimplifation, better use GIVENS ROTATION 
    Jacobian = opD(v);
    u = NuclearNorm(Jacobian); %||J(u)||_nuclear
    
return


function val = NuclearNorm(y)
	s = diff(sum(y.^2,3),1,4);
	theta = atan2(2*dot(y(:,:,:,1),y(:,:,:,2),3),-s)/2;
	c = cos(theta);
	s = sin(theta);
    
    val = sum( sum(sqrt( sum( ( y(:,:,:,1).*c + y(:,:,:,2).*s ).^2, 3) ), 2), 1 )+...
	      sum( sum(sqrt( sum( ( y(:,:,:,2).*c - y(:,:,:,1).*s ).^2, 3) ), 2), 1 );
return

function x = prox_g(y, lambda)

    % GIVENS ROTATION (instead of SVD decomposition)
	s = diff( sum(y.^2,3), 1, 4 );
    
	theta = atan2( 2*dot(y(:,:,:,1), y(:,:,:,2),3), -s )/2;

	c = cos(theta);
	s = sin(theta);

    % S_matrix = [c -s; 
    %             s  c];
	x = cat( 4, y(:,:,:,1).*c + y(:,:,:,2).*s, ...
		       -y(:,:,:,1).*s + y(:,:,:,2).*c ); % x diagonalization 
    
	% PROXIMAL OPERATOR l2 ball
	tmp = max( sqrt(sum(x.^2,3)), lambda );

	tmp =  x .* ( (tmp-lambda)./tmp );

    % S_matrix_T = [c  s;
    %               -s c]
    %     xx = tmp*S_matrix_T 

	x = cat(4, tmp(:,:,:,1).*c - tmp(:,:,:,2).*s, ...
		       tmp(:,:,:,1).*s + tmp(:,:,:,2).*c );
return