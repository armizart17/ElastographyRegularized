function [SWS_CWT,wt_max] = process_CWT(sono_filt_mov,Properties,rd)
%   Function that return a SWS video using the Continuous Wavelet Transform
%   Uses zero padding.
%   Author: Sebastian Merino
%
%   Inputs:
%       sono_filt_mov       Filtered sonoelasticity video
%       Properties          Sonoelasticity video properties, VibFreq and pitch
%       rd                  SWS range
%
%   Outputs:
%       SWS_CWT             SWS video
%       wt_max              Max CWT coefficients

f_lim = 2*Properties.VibFreq./[rd(2),rd(1)]; % SWS range
[f_max,wt_max] = freq_CWT(sono_filt_mov,1/Properties.pitch,f_lim);
SWS_CWT = 2*Properties.VibFreq ./ f_max;