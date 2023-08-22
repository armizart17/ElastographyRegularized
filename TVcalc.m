function [TV] = TVcalc(B,M,N,mask,isotropic)
% Returns the Total Variation (isotropic)
% Inputs: 
%       B               vector containing an image
%       M,N             image size
%       mask            binary mask of analyzed region
%       isotropic       true if the TV is isotropic, false if anisotropic
%  
% Outputs:
%       TV              number
%
% Author: Andres Leonel Coila
% Modified by Sebastian Merino


mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
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