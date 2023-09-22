function [sono_filt_mov,sono_filt,filt_mov]=process_sono_data(sono,Properties,move,rd)

% Median filter
for t = 1:size(sono,3)
    sono(:,:,t) = medfilt2(sono(:,:,t), [3,3], 'symmetric');
end

% Depth normalization
sono_filt = zeros(size(sono));
for id_z = 1:size(sono,1)
    slice = sono(id_z,:,:);
    sono_filt(id_z,:,:) = (slice - mean(slice(:)))./std(slice(:));
end
%sono_filt = normalize(sono,2);

% 2D directional filter
[sono_filt_mov,filt_mov] = moving_filter2(sono_filt,Properties.FrameRate,Properties.pitch,Properties.VibFreq,Properties.VibFreqOffset,rd,move);