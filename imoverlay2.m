function [hF,hB] = imoverlay2(B,SWS,climB,climSWS,alpha,x,z,ROI)
% IMOVERLAY(B,F) displays the image SWS transparently over the image B.
%   alpha:  transparency
%   x:      lateral coordinate in mm
%   z:      depth in mm
%   ROI:    Region of interest
B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);

hB = imagesc(x,z,B);%axis image on;
axis equal
xlim([x(1) x(end)])
ylim([z(1) z(end)])
xlabel('\bf x [mm]')
ylabel('\bf z [mm]')
title('SWS image'); colormap(gray)

h2 = colorbar; ylabel(h2,'m/s'); colormap turbo;
hold on;
hF = imagesc(SWS,climSWS);
% If images are different sizes, map the front image to back coordinates
set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))
% Make the foreground image transparent
alphadata = alpha.*(ROI);
set(hF,'AlphaData',alphadata);
hold off

%set(gcf,'Visible','on');
