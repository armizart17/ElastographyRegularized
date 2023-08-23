function [x,ite_cgs] = cgs2(A,b,tol,maxit,varargin)
% Solves the system of linear equations A*x = b for x using the 
% CONJUGATE GRADIENTS SQUARED METHOD
% Inputs: 
%       A,b         inputs
%       tol         tolerance for error norm
%       maxit       maximum iterations
%       varargin    initial guess for x
%  
% Outputs:
%       x           column vector containing the answer
%       ite_cgs     number of iterations

if length(varargin) == 1
    x = varargin{1};
else
    x= zeros(size(A,2),1);
end
r = A*x - b;
p = -r;
ite_cgs = 0;

while norm(r,2) > tol && ite_cgs < maxit
    alpha = (r'*r)/(p'*(A*p));
    x = x + alpha*p;
    rn = r + alpha*(A*p);
    beta = (rn'*rn)/(r'*r);
    p = -rn + beta*p;
    r = rn;
    ite_cgs = ite_cgs + 1;
    %fprintf("CGS iteration %d\n",ite_cgs);
end

end