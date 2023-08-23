function [TV] = TVcalc(B,M,N,mask,isotropic,colMajor)
% Returns the Total Variation
% Inputs: 
%       B               vector containing an image
%       M,N             image size
%       mask            binary mask of analyzed region
%       isotropic       true if the TV is isotropic, false if anisotropic
%       colMajor        true if the matrix is stored col-major, false if
%                       stored by rows
%  
% Outputs:
%       TV              number
%
% Author: Andres Leonel Coila
% Modified by Sebastian Merino

mask(isnan(mask)) = 0;
mask = mask(:);

if colMajor
    X = reshape(B,M,N);
else
    X = reshape(B,N,M)';
end

Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

if (isotropic)
    P = Dh.^2 + Dv.^2;
    P = sqrt(P);
    TV = norm(P(:).*mask,1);
else
    P = abs(Dh) + abs(Dv);
    TV = norm(P(:).*mask,1);
end

end