function width = fwhm(x,y)

% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'


y = y / max(y);
dx = x(2)-x(1);
N = length(y);

[~,imax]=max(y);
if imax == 1 || imax == N
    disp('Error');
    width = 0;
    return;
end

for iLeft = imax:-1:1
    if y(iLeft) < 0.5
        xLeft = x(iLeft) + dx*(0.5-y(iLeft))/(y(iLeft+1) - y(iLeft));
        break;
    end
    if (iLeft == 1)
        disp('Error');
    end
end

for iRight = imax:N
    if y(iRight) < 0.5
        xRight = x(iRight) - dx*(0.5-y(iRight))/(y(iRight) - y(iRight-1));
        break;
    end
    if (iRight == N)
        disp('Error');
    end
end
%fprintf('xLeft: %.2f, wRight: %.2f\n',xLeft,xRight);
width = xRight - xLeft;
