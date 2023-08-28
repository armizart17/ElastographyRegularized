% Regularized Wavelength Velocity Estimator (R-WAVE) using iterative
% Tikhonov regularization for the inversion problem of each depth
% Inputs:
%
%   sono_filt:      normalized sonoelasticity data, 3D matrix
%
%   Prop.
%       Freq:      vibration freq for CrW generation           [Hz] 
%       Pitch:     spacing between transducer elements         [m]
%
%   RW_Param.
%       dual:       TRUE = valleys included                    [boolean]
%       k:          norm order                                 [ ]
%       beta:       offset to ensure convex solution           [ ]
%       tolerance:  stopping criteria                          [ ]
%       operator:   I = Identity, G = Gradient                 [ ]
%       alpha:      regularization coefficient                 [ ]
%   
% Outputs:
%
%   WA:             Shear Wave speed image                      % [m/s]
%
% Author:           Eduardo A. Gonzalez

function WA = rWave2(sono_filt,Prop,RW_Param)
WA = zeros(size(sono_filt,1),size(sono_filt,2));

for depth = 1:size(sono_filt,1)
    clear B;
    A = sparse(1,size(sono_filt,2));
    n = 1;
    for frame = 1:size(sono_filt,3)
        % Finding indices for each peak in each line (1x128)
        sonoLine = squeeze(sono_filt(depth,:,frame));
        [iPeaks] = peakfinder(sonoLine/max(sonoLine(:)),0);
        
        % Ignoring if peaks are in the start, the end, or there are none
        if size(iPeaks,2)<2, continue; end
        if(iPeaks(1)==1), iPeaks=iPeaks(2:end); end          
        if(iPeaks(end)==size(sono_filt,2)), iPeaks=iPeaks(1:end-1); end
        if size(iPeaks,2)<2, continue; end
    
        % SWS for each wavelength
        lamdaSamples = diff(iPeaks);
        swsLamda = lamdaSamples*2* Prop.pitch * Prop.VibFreq;
    
        % Extrapolation at the start
        A(n,1:iPeaks(1)) = 1/(iPeaks(1));
        B(n) = swsLamda(1);
        n = n+1;
    
        % Estimation of weighted coefficients
        for i=1:size(swsLamda,2) % For each wavelength
            A(n,iPeaks(i)+1:iPeaks(i+1))= 1/lamdaSamples(i);
            B(n) = swsLamda(i); 
            n = n+1;
        end
        
        % Extrapolation at the end
        A(n,(iPeaks(end)+1):end) = 1/(size(A,2)-(iPeaks(end)+1)+1); 
        B(n) = swsLamda(end); 
        n = n+1;
    end
    
    if RW_Param.dual
        % Find valleys too
        for frame=1:size(sono_filt,3)
            sonoLine=squeeze(-1*sono_filt(depth,:,frame));
    
            % Repeating the same procedure
            [iPeaks] = peakfinder(sonoLine/max(sonoLine(:)),0);
            if size(iPeaks,2)<2, continue; end
            if(iPeaks(1)==1), iPeaks=iPeaks(2:end); end
            if(iPeaks(end)==size(sono_filt,2)), iPeaks=iPeaks(1:end-1); end
            if size(iPeaks,2)<2, continue; end
            lamdaSamples = diff(iPeaks);
            swsLamda = lamdaSamples*2* Prop.pitch * Prop.VibFreq;
            A(n,1:iPeaks(1)) = 1/(iPeaks(1));
            B(n) = swsLamda(1);
            n = n+1;
            for i=1:size(swsLamda,2) % For each wavelength
                A(n,iPeaks(i)+1:iPeaks(i+1))= 1/lamdaSamples(i);
                B(n) = swsLamda(i);
                n = n+1;
            end
            A(n,(iPeaks(end)+1):end) = 1/(size(A,2)-(iPeaks(end)+1)+1);
            B(n) = swsLamda(end);
            n = n+1;
        end
    end
    WA(depth,:)=itreg(A,B,RW_Param);
end

end
%% Utility functions
function x = itreg(A,B,RW_Param)
% Choice of regularization function (Tikhonov matrix)
switch (RW_Param.operator)
    case 'I'        % Penalizing greater SWS values
        L=speye(size(A,2));
    case 'G'        % Penalizing greater gradients
        L=speye(size(A,2))-circshift(speye(size(A,2)),[0 1]);
        L(end,:)=0;
end

B = B';
x0 = zeros(size(A,2),1);
err = 1;
while err > RW_Param.tolerance
    Lx = L*x0;
    W = diag( RW_Param.k/2*( abs(Lx.^2+RW_Param.beta).^(RW_Param.k/2 - 1) ) );
    x = ((A'*A+RW_Param.alpha^2 *L'*W*L)\A') *B;
    err = norm(x-x0)^2/norm(x)^2; 
    x0 = x;
end

end