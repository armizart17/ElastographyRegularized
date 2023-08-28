function [f_max,wt_max] = freq_CWT(video,Fs,f_lim)
%   Function that returns the maximum local frequency using CWT.
%   Author: Sebastian Merino
%
%   Inputs:
%       video               Filtered video/image
%       Fs                  Sampling frequency
%       f_lim               Frequency range
%
%   Outputs:
%       f_max               Maximum frequencies
%       wt_max              Max CWT coefficients

[Nz,Nx,Nt] = size(video);
f_max = zeros(Nz,Nx,Nt);
wt_max = f_max;

fb = cwtfilterbank('WaveletParameters',[3,4], ...
    'SamplingFrequency',Fs, 'SignalLength',Nx*2, ...
    'Boundary','periodic', 'VoicesPerOctave',48, 'FrequencyLimits',f_lim);

parfor t_id = 1:Nt
    % Zero-padding for one side because the signal is periodic
    fk_mat = zeros(Nz,Nx);
    U = [video(:,:,t_id),zeros(Nz,Nx)];
    for i=1:Nz
        [wt_coef,fk] = wt(fb,U(i,:));
        [m,id_fk] = max( wt_coef(:,1:Nx) ); 
        fk_mat(i,:) = fk(id_fk)';
        wt_max(i,:,t_id) = m;
    end
    f_max(:,:,t_id) = fk_mat; 
end
end
