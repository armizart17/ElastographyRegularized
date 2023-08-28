function [A,B] = generateWeightMatrix(sonoSlice,Properties)
% Function that generates a linear system of the form Ax = b from
% a depth slice of a sonoelastography video
% Inputs:
%   sonoSlice:  2D matrix filtered sonoelasticity data, where each
%               column is a frame of the same depth
%   Properties
%       .pitch      Separation between piezoelectric elements [m]
%       .VibFreq    Vibration frequency
%   
% Outputs:
%   A:          WxM matrix of weights, where W is the number of measured 
%               wavelengths in all the array and M is the lateral samples.
%   B:          Vector of SWS measurements, of size Wx1.
%
% Code from Eduardo A. Gonzalez adapted by Sebastian Merino
[Nx,Nt] = size(sonoSlice);
maxW = floor(Nx/2*Nt);
B = zeros(maxW,1);
A = spalloc(maxW,Nx,Nx*Nt);
n = 1;
for it = 1:Nt
    % Finding indices for each peak in each line (1x128)
    sonoLine = sonoSlice(:,it);
    [iPeaks] = peakfinder(sonoLine/max(sonoLine(:)),0);
    
    % Ignoring if peaks are in the start, the end, or there are none
    if length(iPeaks)<2, continue; end
    if(iPeaks(1)==1), iPeaks=iPeaks(2:end); end          
    if(iPeaks(end)==Nx), iPeaks=iPeaks(1:end-1); end
    if length(iPeaks)<2, continue; end

    % SWS for each wavelength
    lamdaSamples = diff(iPeaks);
    swsLamda = lamdaSamples*2* Properties.pitch * Properties.VibFreq;

    % Extrapolation at the start (nearest neighbour)
    A(n,1:iPeaks(1)) = 1/(iPeaks(1));
    B(n) = swsLamda(1);
    n = n+1;

    % Estimation of weighted coefficients
    for i = 1:length(swsLamda) % For each wavelength
        A(n,iPeaks(i)+1:iPeaks(i+1)) = 1/lamdaSamples(i);
        B(n) = swsLamda(i); 
        n = n+1;
    end
    
    % Extrapolation at the end
    A(n,(iPeaks(end)+1):Nx) = 1/(Nx-(iPeaks(end)+1)+1); 
    B(n) = swsLamda(end); 
    n = n+1;
end
A(n:end,:) = [];
B(n:end,:) = [];

end