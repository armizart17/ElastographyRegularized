% TV Andres Leonel Coila
function [TV] = TVcalc(B,M,N,mask)
% Returns the Total Variation (isotropic) 
% Inputs: 
%       B               vector containing an image
%       M,N             image size
%       mask            binary mask of analyzed region
%  
% Outputs:
%       TV              number

mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

P = Dh.^2 + Dv.^2;
P = sqrt(P);
TV = norm(P(:).*mask,1);
end

% TV Andres Leonel Coila - ANISOTROPIC
function [TV] = TVcalc2(B,M,N,mask)
% Returns the Total Variation (anisotropic) 
% Inputs: 
%       B               vector containing an image
%       M,N             image size
%       mask            binary mask of analyzed region
%  
% Outputs:
%       TV              number
mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

P = abs(Dh) + abs(Dv); % I don't get this
%P = sqrt(P);
TV = norm(P(:).*mask,1);

end